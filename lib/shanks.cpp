#include <optional>
#include <stdexcept> // For std::runtime_error
#include <limits>    // For std::numeric_limits

#include "gmpxx.h"

#include "factors.hpp"
#include "shanks.hpp"
#include "util.hpp"

namespace primetools {

std::optional<PrimeFactors<mpz_class>>
SQUFOF(
    const mpz_class& N,
    const size_t Max
)
{
    if (N < 2) {
        return std::nullopt;
    }

    mpz_class P0 = sqrt(N);
    mpz_class P_prev = P0;          // P_{i-1} in the forward sequence, starts as P_0
    mpz_class Q_curr = N - P0 * P0; // Q_i in the forward sequence, starts as Q_1

    if (Q_curr == 0) { // N = P0^2, should have been caught by mpz_perfect_square_p. Defensive.
        return PrimeFactors<mpz_class>::FromPair(P0, P0);
    }

    for (size_t fwd_iter_count = 0; fwd_iter_count < Max; ++fwd_iter_count) {
        // Current state for forward search: P_prev is P_{i-1}, Q_curr is Q_i.

        // Calculate P_i = floor((P0 + P_{i-1}) / Q_i) * Q_i - P_{i-1}
        mpz_class P_curr = primetools::floor_div(P0 + P_prev, Q_curr) * Q_curr - P_prev;

        // Check if Q_i (Q_curr) is a positive perfect square
        if (Q_curr > 0 && mpz_perfect_square_p(Q_curr.get_mpz_t())) {
            mpz_class S = sqrt(Q_curr); // S = sqrt(Q_i)

            // Start Reduction Phase
            mpz_class P_prime_prev = P_curr; // P'_0 = P_i (from forward pass)
            mpz_class Q_prime_curr_red;      // Q'_j in reduction sequence

            if (S == 0) { 
                // This case should ideally not be hit if Q_curr > 0.
                // If it occurs, this S cannot be used; continue forward search.
            } else {
                // Q'_1 = (N - (P'_0)^2) / S
                Q_prime_curr_red = primetools::floor_div(N - P_prime_prev * P_prime_prev, S);

                for (size_t red_iter_count = 0; red_iter_count < Max; ++red_iter_count) {
                    if (Q_prime_curr_red == 0) { // Reduction sequence cannot continue
                        break;
                    }

                    // P'_j = floor((P0 + P'_{j-1}) / Q'_j) * Q'_j - P'_{j-1}
                    mpz_class P_prime_curr_val = primetools::floor_div(P0 + P_prime_prev, Q_prime_curr_red) * Q_prime_curr_red - P_prime_prev;

                    if (P_prime_curr_val == P_prime_prev) { // Termination condition for reduction
                        mpz_class factor = gcd(P_prime_curr_val, N);
                        if (factor > 1 && factor < N) {
                            return PrimeFactors<mpz_class>::FromPair(factor, N / factor); // Factor found!
                        }
                        // Trivial factor (1 or N), this S failed. Break reduction.
                        break; 
                    }
                    
                    // Q'_{j+1} = (N - (P'_j)^2) / Q'_j
                    mpz_class Q_prime_next_red = primetools::floor_div(N - P_prime_curr_val * P_prime_curr_val, Q_prime_curr_red);

                    P_prime_prev = P_prime_curr_val;
                    Q_prime_curr_red = Q_prime_next_red;
                } // End reduction loop
            }
            // If reduction loop finishes without returning a factor, this S failed.
            // The forward search loop will continue to the next iteration.
        }

        // Advance forward search: Calculate Q_{i+1} using P_i and Q_i
        // Q_{i+1} = (N - P_i^2) / Q_i
        if (Q_curr == 0) { // Should not happen if Q_curr was used for square check and > 0
            break; // Safety break from forward loop
        }
        mpz_class Q_next = primetools::floor_div(N - P_curr * P_curr, Q_curr);

        // Update P_{i-1} to P_i, and Q_i to Q_{i+1} for the next forward iteration
        P_prev = P_curr;
        Q_curr = Q_next;

        if (Q_curr == 0) { // Forward sequence terminated
            break;
        }
    }

    return std::nullopt; // No factor found within Max
}

} // namespace primetools