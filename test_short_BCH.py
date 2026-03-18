import bchlib
import zlib
import struct
import random
import numpy as np
import matplotlib.pyplot as plt

# --- HARDWARE PARAMETERS ---
DATA_SIZE = 512  # Data packet size (Payload)
OOB_SIZE = 69    # Out-Of-Band physical space
CRC_SIZE = 4     # 32-bit Checksum
a = 2
b = 8
c = 1e5

# FUNCTIONS
def crc32_le(data: bytes) -> bytes:
    """Calculates the 32-bit CRC of the data."""
    return struct.pack("<I", zlib.crc32(data) & 0xFFFFFFFF)

def encode_sector(data: bytes):
    """Write Phase: Calculates CRC, then protects [Data + CRC] with BCH."""
    crc = crc32_le(data)
    protected_payload = data + crc           # Combine data and digital fingerprint
    ecc = bch.encode(protected_payload)      # BCH calculates parity on BOTH
    
    # Build the physical OOB (ECC + CRC + Padding)
    oob = ecc + crc 
    oob += b"\xFF" * (OOB_SIZE - len(oob)) # Padding with 0xFF to reach the OOB size, 
    # the nomenclature is "FF" because in NAND flash, unprogrammed bits are represented by 1s (0xFF per byte)
    return data, oob

def decode_sector(data: bytes, oob: bytes):
    """Read Phase: BCH corrects [Data + CRC], then the CRC validates everything."""
    # 1. Extraction from OOB
    ecc = oob[:bch.ecc_bytes]
    stored_crc = oob[bch.ecc_bytes : bch.ecc_bytes + CRC_SIZE]
    
    # 2. Prepare mutable buffers for BCH correction
    payload_buf = bytearray(data + stored_crc)
    ecc_buf = bytearray(ecc)
    
    # 3. The Surgeon (BCH) detects and corrects errors
    n_errors = bch.decode(payload_buf, ecc_buf)
    if n_errors > 0:
        bch.correct(payload_buf, ecc_buf)
        
    # 4. Separate the repaired data from the repaired CRC
    corrected_data = bytes(payload_buf[:-CRC_SIZE])
    corrected_crc = bytes(payload_buf[-CRC_SIZE:])
    
    # 5. The Inspector (CRC) verifies the surgeon's work
    if n_errors < 0 or crc32_le(corrected_data) != corrected_crc:
        return corrected_data, "UNCORRECTABLE", n_errors
        
    status = "CORRECTED_OK" if n_errors > 0 else "OK"
    return corrected_data, status, n_errors





results = {}

for T in range(a, b,2):
    bch = bchlib.BCH(T, m=14) 
    bit_flips = 1 
    
    for n in range(int(c)):
        if n % (c / 20) == 0 and n > 0:
            bit_flips += 1
            
        original_data = bytes(random.choices(range(256), k=DATA_SIZE))
        
        # 1. Encode
        data_out, oob_out = encode_sector(original_data)
        
        # 2. Corrupt (Dati + OOB)
        full_packet = bytearray(data_out + oob_out)
        for _ in range(bit_flips):
            # pos = random.randint(0, len(full_packet) - 1)
            # full_packet[pos] ^= (1 << random.randint(0, 7)) 
            pos = random.randint(0, len(full_packet) - 1)
            full_packet[pos] ^= 1 
            
        # 3. Decode 
        corrupted_data_part = bytes(full_packet[:DATA_SIZE])
        corrupted_oob_part = bytes(full_packet[DATA_SIZE:])
        
        recovered_data, status, errs = decode_sector(corrupted_data_part, corrupted_oob_part)
        
        # 4. Storage
        key = (T, bit_flips)
        if key not in results: results[key] = []
        results[key].append(1 if status != "UNCORRECTABLE" else 0)

import matplotlib.pyplot as plt
import seaborn as sns # Optional, but makes heatmaps much easier
import numpy as np

# ---------------------------------------------------------
# 1. PROCESS THE DATA
# Calculate the mean success rate (%) for each (T, bit_flips) pair
# ---------------------------------------------------------
success_rates = {}
for key, outcomes in results.items():
    success_rates[key] = np.mean(outcomes) * 100

# Extract unique T values and bit-flips to setup the axes
T_values = sorted(list(set([k[0] for k in success_rates.keys()])))
bit_flip_values = sorted(list(set([k[1] for k in success_rates.keys()])))

# ---------------------------------------------------------
# PLOT 1: Line Chart (Success Rate vs. Bit-flips)
# ---------------------------------------------------------
plt.figure(figsize=(10, 6))
for T in T_values:
    # Get the rates for this specific T, ordered by bit_flips
    rates = [success_rates.get((T, bf), 0) for bf in bit_flip_values]
    plt.plot(bit_flip_values, rates, marker='o', linewidth=2, label=f'T={T}')

plt.axhline(50, color='red', linestyle='--', alpha=0.5, label='50% Threshold')
plt.title('BCH Correction Success Rate vs. Injected Bit-flips', fontsize=14, fontweight='bold')
plt.xlabel('Number of Injected Bit-flips', fontsize=12)
plt.ylabel('Success Rate (%)', fontsize=12)
plt.xticks(bit_flip_values)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()
plt.tight_layout()
plt.show()

# # ---------------------------------------------------------
# # PLOT 2: Grouped Bar Chart
# # ---------------------------------------------------------
# plt.figure(figsize=(12, 6))
# bar_width = 0.8 / len(T_values)
# x = np.arange(len(bit_flip_values))

# for idx, T in enumerate(T_values):
#     rates = [success_rates.get((T, bf), 0) for bf in bit_flip_values]
#     offset = (idx - len(T_values)/2) * bar_width + bar_width/2
#     plt.bar(x + offset, rates, width=bar_width, label=f'T={T}')

# plt.title('Success Rate Histogram by Bit-flips', fontsize=14, fontweight='bold')
# plt.xlabel('Number of Injected Bit-flips', fontsize=12)
# plt.ylabel('Success Rate (%)', fontsize=12)
# plt.xticks(x, bit_flip_values)
# plt.legend(title="Correction Capacity")
# plt.grid(axis='y', linestyle='--', alpha=0.7)
# plt.tight_layout()
# plt.show()

# ---------------------------------------------------------
# PLOT 3: Heatmap (T vs Bit-flips)
# ---------------------------------------------------------
# Create a 2D matrix for the heatmap
rate_matrix = np.zeros((len(T_values), len(bit_flip_values)))

for r, T in enumerate(T_values):
    for c, bf in enumerate(bit_flip_values):
        rate_matrix[r, c] = success_rates.get((T, bf), 0)

plt.figure(figsize=(10, 5))
plt.imshow(rate_matrix, aspect='auto', cmap='RdYlGn', origin='lower', vmin=0, vmax=100)

# Add text annotations inside the heatmap squares
for r in range(len(T_values)):
    for c in range(len(bit_flip_values)):
        val = rate_matrix[r, c]
        color = "white" if val < 30 or val > 70 else "black" # Contrast adjustment
        plt.text(c, r, f"{val:.1f}%", ha='center', va='center', color=color, fontsize=9)

plt.colorbar(label='Success Rate (%)')
plt.title('Heatmap: Success Rate Matrix', fontsize=14, fontweight='bold')
plt.ylabel('Capacity (T)', fontsize=12)
plt.xlabel('Number of Injected Bit-flips', fontsize=12)

plt.yticks(range(len(T_values)), T_values)
plt.xticks(range(len(bit_flip_values)), bit_flip_values)
plt.tight_layout()
plt.show()

