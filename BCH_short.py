import bchlib
import zlib
import struct
import random

# --- HARDWARE PARAMETERS ---
DATA_SIZE = 512  # Data packet size (Payload)
OOB_SIZE = 69    # Out-Of-Band physical space
CRC_SIZE = 4     # 32-bit Checksum
T = 10           # Correction capability (SLC max multi-bit upset)

# Initialize the mathematical universe GF(2^14)
bch = bchlib.BCH(T, m=14) 

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

# ==========================================
# QUICK TEST (Injection of 8 errors)
# ==========================================
if __name__ == "__main__": # This block will run when the script is executed directly
    # Generate 1162 random bytes
    original_data = bytes(random.choices(range(256), k=DATA_SIZE))
    
    # Encode
    data_out, oob_out = encode_sector(original_data)
    
    # Corrupt the data by inserting 8 random bit-flips
    corrupted_data = bytearray(data_out)
    for _ in range(8):
        pos = random.randint(0, len(corrupted_data) - 1)
        corrupted_data[pos] ^= 1  # XOR to invert the bit
        
    # Decode
    recovered_data, status, errs = decode_sector(bytes(corrupted_data), oob_out)
    
    print(f"Status: {status}")
    print(f"Errors found and corrected: {errs}")
    print(f"Is the data perfect? {recovered_data == original_data}")