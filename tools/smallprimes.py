import argparse
import itertools
import multiprocessing
import sys
from Cryptodome.Util.number import isPrime

def is_prime(n):
    """Check if a number _n_ is prime."""
    if n < 2 or (n > 2 and n % 2 == 0):
        return False
    if n in (2, 3,):
        return True
    if n % 3 == 0:
        return False
    for i in range(5, int(n**0.5) + 1, 6):
        if n % i == 0 or n % (i + 2) == 0:
            return False
    return True

def get_primes_in_range(start):
    """Get a list of prime numbers in the range [start, end)."""
    primes = [2,] if start == 0 else []
    end = start + 100_000  # Process in chunks of 100,000
    start = start + 1 if start % 2 == 0 else start
    for num in range(start, end, 2):
        if is_prime(num):
            # print(f"Found prime: {num}")
            primes.append(num)
    return primes

def encode_vlq(value:int):
    """Encode an integer value using Variable Length Quantity (VLQ) encoding."""
    assert value >= 0, "Value must be non-negative"
    while value > 0x7F:
        yield ((value & 0x7F) | 0x80)
        value >>= 7
    yield (value & 0x7F)

def generate_primes_multithreaded(count:int, limit:int = None):
    """Generate a list of the first 'count' prime numbers using multiple threads."""
    found = 0

    # Use a multiprocessing pool to generate blocks of primes
    # in chunks of size 100,000 until the length of primes reaches 'count'.
    chunk_size = 100_000
    start = 0

    last_prime = 0
    limit_reached = False

    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        while found < count and (limit is not None and not limit_reached):
            end = start + chunk_size * multiprocessing.cpu_count()
            results = pool.map(get_primes_in_range, range(start, end, chunk_size))
            primes = list(itertools.chain.from_iterable(results))
            primes.sort()
            assert len(primes) == len(set(primes)), "Primes list contains duplicates"
            for prime in primes:
                if limit is not None and prime >= limit:
                    if not limit_reached:
                        sys.stderr.write(f"\nReached limit of {limit:,d} primes.\n")
                    limit_reached = True
                    break
                if found < count:
                    assert prime > last_prime, f"Expected {prime} to be greater than {last_prime}"
                    assert prime - last_prime > 0, f"Expected {prime} - {last_prime} > 0"
                    value = prime - last_prime
                    yield from encode_vlq(value)
                    last_prime = prime
                    found += 1
                else:
                    break
            sys.stderr.write(f"\rGenerated {found:,d} primes. Last prime: {last_prime:,d}")
            sys.stderr.flush()
            start = end

def generate_primes_linear(limit):
    """Generator that yields the first 'limit' prime numbers."""
    last_prime = 2
    yield last_prime  # First prime
    num = 3
    count = 1
    while count < limit:
        if is_prime(num):
            assert num > last_prime, f"Expected {num} to be greater than {last_prime}"
            assert num - last_prime > 0
            value = num - last_prime
            yield from encode_vlq(value)
            last_prime = num
            count += 1
        num += 2

def output_primes(prime_count: int, file_name: str, limit: int = None):
    """Output the first 'prime_count' primes in a binary file."""
    generated = 0
    
    with open(file_name, "wb") as f:
        for prime in generate_primes_multithreaded(prime_count, limit):
            f.write(prime.to_bytes(1, 'big'))
            generated += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a list of prime numbers.")
    parser.add_argument("count", type=int, help="Number of prime numbers to generate")
    parser.add_argument("file", type=str, default="primes.bin", help="Output file name")
    parser.add_argument("--limit", type=int, help="Limit the maximum prime number to generate")
    
    args = parser.parse_args()
    
    if args.count <= 0:
        print("Count must be a positive integer.")
        sys.exit(1)
    
    output_primes(args.count, args.file, args.limit)
