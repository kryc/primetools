import itertools

HEAD = """#include <array>

static const std::array<unsigned char, {:d}> g_FirstPrimes = {{
"""

FOOT = """
};"""

def is_prime(n):
    if n < 2 or (n > 2 and n % 2 == 0):
        return False
    if n == 2 or n == 3:
        return True
    if n % 3 == 0:
        return False
    for i in range(5, int(n**0.5) + 1, 6):
        if n % i == 0 or n % (i + 2) == 0:
            return False
    return True

def generate_primes(limit):
    """Generator that yields the first 'limit' prime numbers."""
    last_prime = 2
    yield last_prime  # First prime
    num = 3
    count = 1
    while count < limit:
        if is_prime(num):
            assert num > last_prime, f"Expected {num} to be greater than {last_prime}"
            assert num - last_prime > 0
            assert num - last_prime < 256
            yield num - last_prime
            last_prime = num
            count += 1
        num += 2

def generate_primes_list(limit, block_size=32):
    """Generate a list of the first 'limit' prime numbers in blocks of 'block_size'."""
    primes = []
    for prime in itertools.islice(generate_primes(limit), limit):
        primes.append(prime)
        if len(primes) % block_size == 0:
            yield primes
            primes = []
    if primes:
        yield primes

def output_primes(prime_count):
    """Output the first 'prime_count' primes in a C++ array format."""
    assert prime_count > 0
    assert prime_count % 32 == 0
    print(HEAD.format(prime_count))
    for block in generate_primes_list(prime_count):
        line = ", ".join(f"0x{p:02x}" for p in block)
        print(f"    {line},")
    print(FOOT)

if __name__ == "__main__":
    primes = output_primes(3_200_000)
