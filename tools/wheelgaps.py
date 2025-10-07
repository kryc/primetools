import argparse
from Cryptodome.Util import number
from math import gcd
import os
import struct

def calculate_wheel_gaps(max_prime):
    # Get all primes up to and uncluding max_prime
    primes = [p for p in range(2, max_prime + 1) if number.isPrime(p)]
    W = 1
    for p in primes:
        W *= p

    residues = [r for r in range(1, W) if gcd(r, W) == 1]
    gaps = [(residues[(i+1) % len(residues)] - residues[i]) % W
            for i in range(len(residues))]

    # assert all gaps are even
    assert all(g % 2 == 0 for g in gaps)

    # Assert the modulus is even
    assert W % 2 == 0

    # return modulus and gaps
    return W, gaps

def main(out_file: str, max_prime: int, pack: bool = True, use_128: bool = False):
    W, gaps = calculate_wheel_gaps(max_prime)

    # Calculate the number of bits needed to store the largest gap
    largest_gap = max(gaps)
    bits = largest_gap.bit_length()

    storage_type = "uint64_t" if not use_128 else "__uint128_t"
    storage_bits = 64 if not use_128 else 128

    if pack:
        # Shift all gaps right by 1
        gaps = [g >> 1 for g in gaps]
        bits -= 1  # one less bit needed
        storage_bits -= 1  # one less bit available

    bit_mask = (1 << bits) - 1
    gaps_per_word = storage_bits // bits
    word_count = (len(gaps) + gaps_per_word - 1) // gaps_per_word  # Ceiling division

    # 2. Verify N-bit fit
    assert max(gaps) <= bit_mask

    # 3. Emit C‐style packed array
    with open(out_file, "wb") as f:
        print(f"// Wheel modulus = {W}, total gaps = {len(gaps)}, bits_required = {bits}, gaps_per_word = {gaps_per_word}, word_count = {word_count}")
        print(f"static const {storage_type} WHEEL{W}GAPS[] = " + "{")
        print(f"#embed \"{out_file}\"")
        for i in range(0, len(gaps), gaps_per_word):
            block = gaps[i:i+gaps_per_word]
            word = 0
            for j, g in enumerate(block):
                shift = j * bits  # pack each gap
                assert (g & bit_mask) == g
                word |= g << shift
            # Final shift left by 1
            if pack:
                word <<= 1
            if use_128:
                # Split into two 64-bit parts for const assignnemt
                upper = (word >> 64) & 0xFFFFFFFFFFFFFFFF
                lower = word & 0xFFFFFFFFFFFFFFFF
                # print(f"    ((__uint128_t)0x{upper:016x}ULL << 64) | 0x{lower:016x}ULL,")
                f.write(struct.pack("<QQ", lower, upper))
            else:
                # print(f"    0x{word:016x}ULL,")
                f.write(struct.pack("<Q", word))
        print("};")
        print()

    # Validate the file size
    expected_size = word_count * (16 if use_128 else 8)
    actual_size = os.path.getsize(out_file)
    assert expected_size == actual_size, f"Expected file size {expected_size}, got {actual_size}"

if __name__ == "__main__":
    # First primes: 2, 3, 5, 7(30), 11(210), 13(30030), 17(510510), 19(9699690), 23(223092870), 29, 31, 37, 41, 43, 47
    parser = argparse.ArgumentParser(description="Calculate wheel gaps for prime sieving.")
    parser.add_argument("output", type=str, help="Output file name.")
    parser.add_argument("max_prime", type=int, help="Maximum prime to use for wheel calculation.")
    parser.add_argument("--nopack", action="store_true", help="Do not pack gaps, use full size.")
    parser.add_argument("--use128", action="store_true", help="Use 128-bit storage type.")
    args = parser.parse_args()

    main(args.output, args.max_prime, pack=not args.nopack, use_128=args.use128)