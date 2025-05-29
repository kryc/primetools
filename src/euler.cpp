#include "euler.hpp"

const mpz_class
EulerTotient(
    const mpz_class N
)
{
    if (N < 2) {
        return 0;
    }

    mpz_t result;
    mpz_init_set(result, N.get_mpz_t()); // Initialize result to N
    mpz_t p, n, rootN, temp;
    mpz_inits(n, rootN, temp);

    mpz_set(n, N.get_mpz_t());
    mpz_sqrt(rootN, n); // Calculate sqrt(N)

    mpz_init_set_ui(p, 2); // Start with the smallest prime factor

    // Check for factors up to sqrt(N)
    while (mpz_cmp(p, rootN) <= 0) {
        if (mpz_divisible_p(n, p)) {
            // p is A prime factor of N
            while (mpz_divisible_p(n, p)) {
                mpz_divexact(n, n, p); // Divide N by p
            }
            mpz_divexact(temp, result, p); // Calculate result / p
            mpz_sub(result, result, temp); // result = result - (result / p)
        }
        // Move to the next prime factor
        mpz_add_ui(p, p, 1); // Increment p
    }

    // If N is still greater than 1, then it is prime
    if (mpz_cmp_ui(n, 1) > 0) {
        mpz_divexact(temp, result, n); // Calculate result / N
        mpz_sub(result, result, temp); // result = result - (result / N)
    }

    return mpz_class(result);
}

const bool
EulerCriterion(
    const mpz_class A,
    const mpz_class P
)
{
    if (P <= 1 || A < 0 || A >= P) {
        return false; // Invalid input
    }

    mpz_t exponent, result;
    mpz_inits(exponent, result);

    // Calculate A^(P-1)/2 mod P
    mpz_sub_ui(exponent, P.get_mpz_t(), 1); // exponent = P - 1
    mpz_div_ui(exponent, exponent, 2); // exponent = (P - 1) / 2
    mpz_powm(result, A.get_mpz_t(), exponent, P.get_mpz_t());

    return mpz_cmp_ui(result, 1); // If result is 1, then A is A quadratic residue mod P
}