# Primetools
By Kryc

## Introduction

Primetools is a set of utilities for manipulating prime numbers and, more specifically, factorising the product of two primes. The latter is extremely useful for attacking weak RSA keys.

The main use of this project is in the building and solving of capture the flag (CTF) competitions.

## Factorisation Strategies

- Small prime factorisation (against the first 100K primes)
- [Fermat's factorisation method](https://en.wikipedia.org/wiki/Fermat's_factorization_method)
- [Pollard's Rho](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm)
- [Brent-Pollards Rho](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants)

Note that although the original Pollard's Rho is implemented and working, Brent's improvements provide a roughly 24% speed improvement so this version is used as the default.

## Building

Building the project is straightforward on *nix systems.

```bash
mkdir build;
cd build;
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Running

There are a couple of command line switches, but the main one to use would be:

```bash
./primetools factorise <N>
```

Where `N` is the number whose prime factors you wish to find.
