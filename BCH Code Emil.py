import bchlib
import struct
import zlib
import random
from collections import defaultdict
import matplotlib.pyplot as plt

# --------------------------
# Sector & Overhead Definition
# --------------------------
SECTOR_SIZE = 512
OOB_PER_SECTOR = 69
CRC_SIZE = 4  # CRC-32

# ---------------------------
# 1) BCH initialization
# ---------------------------
T = 4

def make_bch(m: int) -> bchlib.BCH:
    bch = bchlib.BCH(T, m=m)
    if bch.ecc_bytes + CRC_SIZE > OOB_PER_SECTOR:
        raise ValueError(
            f"ECC too large: ecc_bytes={bch.ecc_bytes}, "
            f"need <= {OOB_PER_SECTOR - CRC_SIZE} bytes for ECC"
        )
    return bch

bch = make_bch(m=13)
ECC_BYTES = bch.ecc_bytes
print("BCH configured:", {"m": 13, "t": T, "ecc_bytes": ECC_BYTES})

# ---------------------------
# Helpers
# ---------------------------
def crc32_le(data: bytes) -> bytes:
    c = zlib.crc32(data) & 0xFFFFFFFF
    return struct.pack("<I", c)

def pack_oob(ecc: bytes, crc: bytes) -> bytes:
    assert len(crc) == CRC_SIZE
    assert len(ecc) == ECC_BYTES
    oob = ecc + crc
    oob += b"\xFF" * (OOB_PER_SECTOR - len(oob))
    return oob

def unpack_oob(oob: bytes):
    if len(oob) != OOB_PER_SECTOR:
        raise ValueError(f"OOB must be exactly {OOB_PER_SECTOR} bytes, got {len(oob)}")
    ecc = oob[:ECC_BYTES]
    crc = oob[ECC_BYTES:ECC_BYTES + CRC_SIZE]
    return ecc, crc

# ---------------------------
# 2) Sector encode / decode
# ---------------------------
def encode_sector(data512: bytes):
    if len(data512) != SECTOR_SIZE:
        raise ValueError(f"sector must be {SECTOR_SIZE} bytes, got {len(data512)}")
    ecc = bch.encode(data512)
    crc = crc32_le(data512)
    oob = pack_oob(ecc, crc)
    return data512, oob

class DecodeStatus:
    OK = "OK"
    CORRECTED_OK = "CORRECTED_OK"
    UNCORRECTABLE = "UNCORRECTABLE"

def decode_sector(data512: bytes, oob69: bytes):
    if len(data512) != SECTOR_SIZE:
        raise ValueError(f"data must be {SECTOR_SIZE} bytes")
    if len(oob69) != OOB_PER_SECTOR:
        raise ValueError(f"OOB must be {OOB_PER_SECTOR} bytes")

    stored_ecc, stored_crc = unpack_oob(oob69)

    data_buf = bytearray(data512)
    ecc_buf = bytearray(stored_ecc)

    try:
        n_errors = bch.decode(data_buf, ecc_buf)
    except Exception:
        return bytes(data512), DecodeStatus.UNCORRECTABLE, 0

    if n_errors < 0:
        return bytes(data512), DecodeStatus.UNCORRECTABLE, 0

    if n_errors > 0:
        bch.correct(data_buf, ecc_buf)

    corrected = bytes(data_buf)

    calc_crc = crc32_le(corrected)
    if calc_crc != stored_crc:
        return corrected, DecodeStatus.UNCORRECTABLE, n_errors

    status = DecodeStatus.CORRECTED_OK if n_errors > 0 else DecodeStatus.OK
    return corrected, status, n_errors

# ---------------------------
# 3) Fault injection
# ---------------------------
def flip_bits(buf: bytes, bit_positions: list) -> bytes:
    b = bytearray(buf)
    for pos in bit_positions:
        byte_i = pos // 8
        bit_i = pos % 8
        b[byte_i] ^= (1 << bit_i)
    return bytes(b)

def random_flip_bits(buf: bytes, n_flips: int) -> bytes:
    total_bits = len(buf) * 8
    positions = random.sample(range(total_bits), n_flips)
    return flip_bits(buf, positions)

# ---------------------------
# 4) Single trial
# ---------------------------
def run_trial(n_data_flips=0, n_ecc_flips=0, n_crc_flips=0, seed=None):
    if seed is not None:
        random.seed(seed)

    original = bytes([i % 256 for i in range(SECTOR_SIZE)])
    data_enc, oob = encode_sector(original)

    if n_data_flips > 0:
        data_enc = random_flip_bits(data_enc, n_data_flips)

    ecc, crc = unpack_oob(oob)

    if n_ecc_flips > 0:
        ecc = random_flip_bits(ecc, n_ecc_flips)

    if n_crc_flips > 0:
        crc = random_flip_bits(crc, n_crc_flips)

    oob_faulty = pack_oob(ecc, crc)

    corrected, status, n_corr = decode_sector(data_enc, oob_faulty)
    return {
        "status": status,
        "n_corrected": n_corr,
        "data_matches_original": (corrected == original)
    }

# ---------------------------
# 5) Monte Carlo campaign
# ---------------------------
def monte_carlo(flips, trials):
    counts = defaultdict(int)
    corrected_hist = defaultdict(int)
    exact_match = 0

    for i in range(trials):
        r = run_trial(n_data_flips=flips, seed=i)

        counts[r["status"]] += 1
        corrected_hist[r["n_corrected"]] += 1

        if r["data_matches_original"]:
            exact_match += 1

    result = {
        "flips": flips,
        "trials": trials,
        "OK": counts[DecodeStatus.OK],
        "CORRECTED_OK": counts[DecodeStatus.CORRECTED_OK],
        "UNCORRECTABLE": counts[DecodeStatus.UNCORRECTABLE],
        "exact_match": exact_match,
        "failure_rate": counts[DecodeStatus.UNCORRECTABLE] / trials,
        "success_rate": exact_match / trials,
        "corrected_hist": dict(corrected_hist),
    }
    return result

def print_results_table(results):
    print("\n=== Monte Carlo Summary ===")
    print(f"{'Flips':>5} | {'Trials':>8} | {'OK':>8} | {'CORR_OK':>10} | {'UNCORR':>8} | {'Success %':>10} | {'Failure %':>10}")
    print("-" * 80)
    for r in results:
        print(
            f"{r['flips']:>5} | "
            f"{r['trials']:>8} | "
            f"{r['OK']:>8} | "
            f"{r['CORRECTED_OK']:>10} | "
            f"{r['UNCORRECTABLE']:>8} | "
            f"{100*r['success_rate']:>9.4f}% | "
            f"{100*r['failure_rate']:>9.4f}%"
        )

    print("\n=== Correction Count Histograms ===")
    for r in results:
        print(f"\n{r['flips']}-bit flip trials:")
        for k in sorted(r["corrected_hist"].keys()):
            print(f"  n_corrected = {k}: {r['corrected_hist'][k]}")

def plot_results(results):
    flip_labels = [str(r["flips"]) for r in results]
    success_rates = [100 * r["success_rate"] for r in results]
    failure_rates = [100 * r["failure_rate"] for r in results]

    plt.figure(figsize=(10, 6))
    plt.bar(flip_labels, success_rates)
    plt.ylabel("Exact recovery rate (%)")
    plt.xlabel("Injected data bit flips")
    plt.title("BCH+CRC recovery effectiveness")
    plt.ylim(0, 105)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.bar(flip_labels, failure_rates)
    plt.ylabel("Failure rate (%)")
    plt.xlabel("Injected data bit flips")
    plt.title("BCH+CRC failure rate")
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.show()

# ---------------------------
# 6) Main campaign
# ---------------------------
if __name__ == "__main__":
    campaign = [
        {"flips": 1, "trials": 2000000},
        {"flips": 2, "trials": 20000000},
        {"flips": 3, "trials": 10000},
    ]

    results = []
    for case in campaign:
        print(f"Running {case['trials']} trials for {case['flips']}-bit flips...")
        results.append(monte_carlo(case["flips"], case["trials"]))

    print_results_table(results)
    plot_results(results)