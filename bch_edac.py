"""
BCH EDAC simulator — improved with live dashboard, verbose tracing, and 2-bit exhaustive hunt.

Usage:
  python bch_edac.py --demo                        # verbose single-sector walkthrough
  python bch_edac.py --run [--flips N] [--trials N] # live-dashboard Monte Carlo
  python bch_edac.py --hunt2                        # exhaustive 2-bit uncorrectable search
"""

import bchlib
import struct
import zlib
import random
import sys
import time
import json
import argparse
import math
from collections import defaultdict

# ─────────────────────────────────────────────
# ANSI helpers
# ─────────────────────────────────────────────
class C:
    GREEN  = "\033[32m"
    YELLOW = "\033[33m"
    RED    = "\033[31m"
    CYAN   = "\033[36m"
    WHITE  = "\033[97m"
    DIM    = "\033[2m"
    BOLD   = "\033[1m"
    RESET  = "\033[0m"

def ok(s):   return f"{C.GREEN}{s}{C.RESET}"
def warn(s): return f"{C.YELLOW}{s}{C.RESET}"
def err(s):  return f"{C.RED}{s}{C.RESET}"
def dim(s):  return f"{C.DIM}{s}{C.RESET}"
def bold(s): return f"{C.BOLD}{s}{C.RESET}"

def progress_bar(done, total, width=30):
    frac = done / total if total else 0
    filled = int(frac * width)
    bar = "█" * filled + "░" * (width - filled)
    return f"[{bar}] {frac*100:5.1f}%"

# ─────────────────────────────────────────────
# BCH / sector constants
# ─────────────────────────────────────────────
SECTOR_SIZE    = 512
OOB_PER_SECTOR = 69
CRC_SIZE       = 4
T              = 4
BCH_M          = 13

bch = bchlib.BCH(T, m=BCH_M)
ECC_BYTES = bch.ecc_bytes

# ─────────────────────────────────────────────
# Codec helpers
# ─────────────────────────────────────────────
def crc32_le(data: bytes) -> bytes:
    c = zlib.crc32(data) & 0xFFFFFFFF
    return struct.pack("<I", c)

def pack_oob(ecc: bytes, crc: bytes) -> bytes:
    oob = ecc + crc
    oob += b"\xFF" * (OOB_PER_SECTOR - len(oob))
    return oob

def unpack_oob(oob: bytes):
    ecc = oob[:ECC_BYTES]
    crc = oob[ECC_BYTES:ECC_BYTES + CRC_SIZE]
    return ecc, crc

def encode_sector(data: bytes):
    ecc = bch.encode(data)
    crc = crc32_le(data)
    oob = pack_oob(ecc, crc)
    return data, oob

class Status:
    OK           = "OK"
    CORRECTED_OK = "CORRECTED_OK"
    UNCORRECTABLE = "UNCORRECTABLE"

def decode_sector(data: bytes, oob: bytes):
    stored_ecc, stored_crc = unpack_oob(oob)
    data_buf = bytearray(data)
    ecc_buf  = bytearray(stored_ecc)
    try:
        n_errors = bch.decode(data_buf, ecc_buf)
    except Exception:
        return bytes(data), Status.UNCORRECTABLE, 0
    if n_errors < 0:
        return bytes(data), Status.UNCORRECTABLE, 0
    if n_errors > 0:
        bch.correct(data_buf, ecc_buf)
    corrected = bytes(data_buf)
    if crc32_le(corrected) != stored_crc:
        return corrected, Status.UNCORRECTABLE, n_errors
    return corrected, (Status.CORRECTED_OK if n_errors > 0 else Status.OK), n_errors

# ─────────────────────────────────────────────
# Fault injection
# ─────────────────────────────────────────────
def flip_bits(buf: bytes, positions: list) -> bytes:
    b = bytearray(buf)
    for pos in positions:
        b[pos // 8] ^= (1 << (pos % 8))
    return bytes(b)

def random_bit_positions(buf_len: int, n: int) -> list:
    return random.sample(range(buf_len * 8), n)

# ─────────────────────────────────────────────
# MODE 1 — verbose demo
# ─────────────────────────────────────────────
def run_demo(n_flips: int = 2, seed: int = 42):
    random.seed(seed)
    print()
    print(bold("─" * 60))
    print(bold("  BCH EDAC — verbose sector walkthrough"))
    print(bold("─" * 60))
    print(f"  BCH params : m={BCH_M}, t={T} (corrects up to {T} bit flips)")
    print(f"  ECC bytes  : {ECC_BYTES}")
    print(f"  OOB layout : [{ECC_BYTES}B ECC | {CRC_SIZE}B CRC | {OOB_PER_SECTOR-ECC_BYTES-CRC_SIZE}B pad]")
    print(f"  Injecting  : {n_flips} random bit flip(s)")
    print()

    original = bytes([i % 256 for i in range(SECTOR_SIZE)])
    print(ok("  [1] Encoding sector"))
    data_enc, oob = encode_sector(original)
    ecc, crc = unpack_oob(oob)
    print(f"       Data (first 16B) : {original[:16].hex(' ')}")
    print(f"       ECC ({ECC_BYTES}B)       : {ecc.hex(' ')}")
    print(f"       CRC (4B)          : {crc.hex(' ')}")

    print()
    positions = random_bit_positions(SECTOR_SIZE, n_flips)
    positions.sort()
    print(warn(f"  [2] Injecting {n_flips} bit flip(s)"))
    for p in positions:
        byte_i, bit_i = p // 8, p % 8
        before = original[byte_i]
        after  = before ^ (1 << bit_i)
        print(f"       bit {p:5d}  →  byte[{byte_i:3d}] bit{bit_i}  : 0x{before:02X} → 0x{after:02X}")

    corrupted = flip_bits(data_enc, positions)
    diffs = sum(1 for a, b in zip(original, corrupted) if a != b)
    print(f"       Corrupted bytes   : {diffs}")

    print()
    print(f"  [3] Decoding …")
    corrected, status, n_corr = decode_sector(corrupted, oob)

    if status == Status.OK:
        print(ok(f"       Status : OK (no errors detected)"))
    elif status == Status.CORRECTED_OK:
        print(warn(f"       Status : CORRECTED_OK ({n_corr} error(s) fixed)"))
    else:
        print(err(f"       Status : UNCORRECTABLE"))

    matches = corrected == original
    print(f"       Data matches original : {ok('YES') if matches else err('NO')}")

    if not matches:
        print()
        print(err("  [!] Byte-level diff (first 8 mismatches):"))
        shown = 0
        for i, (o, c) in enumerate(zip(original, corrected)):
            if o != c:
                print(f"       byte[{i:3d}] : original=0x{o:02X}  corrected=0x{c:02X}")
                shown += 1
                if shown >= 8:
                    break

    print()
    print(bold("─" * 60))
    print()

# ─────────────────────────────────────────────
# Monte Carlo trial (no I/O, fast)
# ─────────────────────────────────────────────
def run_trial(n_flips: int, seed: int):
    random.seed(seed)
    original = bytes([i % 256 for i in range(SECTOR_SIZE)])
    data_enc, oob = encode_sector(original)
    positions = random_bit_positions(SECTOR_SIZE, n_flips)
    corrupted = flip_bits(data_enc, positions)
    corrected, status, n_corr = decode_sector(corrupted, oob)
    return status, corrected == original, n_corr, positions

# ─────────────────────────────────────────────
# MODE 2 — live-dashboard Monte Carlo run
# ─────────────────────────────────────────────
def run_live(n_flips: int, total_trials: int, log_file: str = "bch_results.json"):
    counts = defaultdict(int)
    failures = []
    t_start = time.time()
    t_last_print = t_start - 1

    print()
    print(bold(f"  BCH Monte Carlo — {n_flips}-bit flips, {total_trials:,} trials"))
    print(dim(f"  Results will be written to {log_file}"))
    print()

    for i in range(total_trials):
        status, match, n_corr, positions = run_trial(n_flips, seed=i)
        counts[status] += 1

        # record any uncorrectable failures (store first 100)
        if status == Status.UNCORRECTABLE and len(failures) < 100:
            failures.append({"seed": i, "bit_positions": sorted(positions)})

        # live dashboard refresh (at most once per second)
        now = time.time()
        if now - t_last_print >= 1.0 or i == total_trials - 1:
            elapsed = now - t_start
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta  = (total_trials - i - 1) / rate if rate > 0 else 0

            n_ok    = counts[Status.OK]
            n_corr_ = counts[Status.CORRECTED_OK]
            n_fail  = counts[Status.UNCORRECTABLE]
            done    = i + 1

            bar = progress_bar(done, total_trials)

            sys.stdout.write("\033[2K\r")  # clear line
            line = (
                f"  {bar}  "
                f"{ok(f'OK:{n_ok:>8,}')}  "
                f"{warn(f'CORR:{n_corr_:>8,}')}  "
                f"{err(f'FAIL:{n_fail:>6,}')}  "
                f"{dim(f'{rate:,.0f} t/s  ETA {eta:.0f}s')}"
            )
            sys.stdout.write(line)
            sys.stdout.flush()
            t_last_print = now

    elapsed = time.time() - t_start
    print()  # newline after live bar
    print()

    # ── summary ──
    n_ok    = counts[Status.OK]
    n_corr_ = counts[Status.CORRECTED_OK]
    n_fail  = counts[Status.UNCORRECTABLE]

    print(bold("  ── Results ─────────────────────────────────"))
    print(f"  Trials       : {total_trials:>12,}")
    print(f"  {ok('OK           :')} {n_ok:>12,}  ({100*n_ok/total_trials:.4f}%)")
    print(f"  {warn('Corrected    :')} {n_corr_:>12,}  ({100*n_corr_/total_trials:.4f}%)")
    print(f"  {err('Uncorrectable:')} {n_fail:>12,}  ({100*n_fail/total_trials:.6f}%)")
    print(f"  Elapsed      : {elapsed:.1f}s  ({total_trials/elapsed:,.0f} trials/s)")
    print()

    if failures:
        print(err(f"  [!] {len(failures)} uncorrectable case(s) logged (first 100 max):"))
        for f in failures[:5]:
            print(f"       seed={f['seed']}  bits={f['bit_positions']}")
        if len(failures) > 5:
            print(dim(f"       … and {len(failures)-5} more in {log_file}"))
    else:
        print(ok("  No uncorrectable cases found."))

    print()

    # ── write JSON log ──
    log = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "n_flips": n_flips,
        "total_trials": total_trials,
        "OK": n_ok,
        "CORRECTED_OK": n_corr_,
        "UNCORRECTABLE": n_fail,
        "failure_rate": n_fail / total_trials,
        "elapsed_s": round(elapsed, 2),
        "uncorrectable_cases": failures,
    }
    with open(log_file, "a") as f:
        f.write(json.dumps(log) + "\n")
    print(dim(f"  Log saved → {log_file}"))
    print()

# ─────────────────────────────────────────────
# MODE 3 — exhaustive 2-bit uncorrectable hunter
# ─────────────────────────────────────────────
def run_hunt2(log_file: str = "bch_hunt2.json"):
    """
    Test every possible pair of bit positions in the 512-byte sector.
    C(4096, 2) = 8,386,560 unique pairs — feasible in a few minutes.
    """
    total_bits = SECTOR_SIZE * 8   # 4096
    total_pairs = total_bits * (total_bits - 1) // 2  # 8,386,560

    print()
    print(bold("  BCH 2-bit exhaustive hunt"))
    print(f"  Testing all C({total_bits}, 2) = {total_pairs:,} unique bit-pair positions")
    print()

    original = bytes([i % 256 for i in range(SECTOR_SIZE)])
    data_enc, oob = encode_sector(original)

    failures = []
    t_start = time.time()
    t_last_print = t_start - 1
    done = 0

    for i in range(total_bits):
        for j in range(i + 1, total_bits):
            corrupted = flip_bits(data_enc, [i, j])
            corrected, status, n_corr = decode_sector(corrupted, oob)

            if status == Status.UNCORRECTABLE:
                failures.append({"bit_a": i, "bit_b": j,
                                  "byte_a": i//8, "byte_b": j//8})

            done += 1
            now = time.time()
            if now - t_last_print >= 1.0 or done == total_pairs:
                elapsed = now - t_start
                rate = done / elapsed if elapsed > 0 else 0
                eta  = (total_pairs - done) / rate if rate > 0 else 0
                bar  = progress_bar(done, total_pairs)
                sys.stdout.write("\033[2K\r")
                sys.stdout.write(
                    f"  {bar}  "
                    f"{err(f'failures: {len(failures):>5}')}"
                    f"  {dim(f'{rate:,.0f} pairs/s  ETA {eta:.0f}s')}"
                )
                sys.stdout.flush()
                t_last_print = now

    elapsed = time.time() - t_start
    print()
    print()

    print(bold("  ── Hunt results ────────────────────────────"))
    print(f"  Pairs tested  : {total_pairs:,}")
    print(f"  Elapsed       : {elapsed:.1f}s  ({total_pairs/elapsed:,.0f} pairs/s)")

    if failures:
        print(err(f"  Uncorrectable : {len(failures):,} pairs ({100*len(failures)/total_pairs:.4f}%)"))
        print()
        print(err("  Sample failures (first 10):"))
        for f in failures[:10]:
            print(f"       bits ({f['bit_a']:4d}, {f['bit_b']:4d})  "
                  f"→  bytes [{f['byte_a']:3d}] [{f['byte_b']:3d}]")
    else:
        print(ok("  All 2-bit patterns correctable — BCH(t=4) is solid for this sector."))

    print()

    log = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "total_pairs": total_pairs,
        "uncorrectable_count": len(failures),
        "failure_rate": len(failures) / total_pairs,
        "elapsed_s": round(elapsed, 2),
        "failures": failures,
    }
    with open(log_file, "a") as f:
        f.write(json.dumps(log) + "\n")
    print(dim(f"  Log saved → {log_file}"))
    print()

# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="BCH EDAC simulator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--demo",   action="store_true", help="verbose single-sector walkthrough")
    group.add_argument("--run",    action="store_true", help="live-dashboard Monte Carlo")
    group.add_argument("--hunt2",  action="store_true", help="exhaustive 2-bit hunt")

    parser.add_argument("--flips",  type=int, default=2,       help="bit flips to inject (--run mode, default 2)")
    parser.add_argument("--trials", type=int, default=5_000_000, help="number of trials (--run mode, default 5M)")
    parser.add_argument("--log",    type=str, default=None,    help="output JSON log filename")

    args = parser.parse_args()

    print()
    print(bold(f"  BCH params: m={BCH_M}  t={T}  ecc_bytes={ECC_BYTES}  "
               f"sector={SECTOR_SIZE}B  OOB={OOB_PER_SECTOR}B"))

    if args.demo:
        run_demo(n_flips=args.flips)

    elif args.run:
        log = args.log or "bch_results.json"
        run_live(n_flips=args.flips, total_trials=args.trials, log_file=log)

    elif args.hunt2:
        log = args.log or "bch_hunt2.json"
        run_hunt2(log_file=log)


if __name__ == "__main__":
    main()
