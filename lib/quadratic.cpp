#include <algorithm>
#include <atomic>
#include <bitset>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <random>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "maths.hpp"
#include "logging.hpp"
#include "primegenerator.hpp"
#include "quadratic.hpp"
#include "tonelli_shanks.hpp"
#include "util.hpp"

namespace primetools {
namespace {

// std::bitset requires compile-time sizes; our QS sizes are runtime.
// We cap them based on current parameter clamping (B <= 200k) and relation slack.
constexpr size_t kMaxCols = 20000;      // sign + primes in factor base
constexpr size_t kMaxRelations = 30000; // relations tracked in elimination

static uint32_t
ModInverseU32(
    uint32_t a,
    uint32_t mod
)
{
    // Extended Euclid for uint32 moduli.
    // Returns x such that (a*x) % mod == 1.
    if (mod == 0) {
        throw std::invalid_argument("ModInverseU32: mod must be non-zero");
    }
    a %= mod;
    if (a == 0) {
        throw std::invalid_argument("ModInverseU32: no inverse for 0");
    }

    int64_t t = 0;
    int64_t newt = 1;
    int64_t r = static_cast<int64_t>(mod);
    int64_t newr = static_cast<int64_t>(a);

    while (newr != 0) {
        const int64_t q = r / newr;
        const int64_t tmp_t = t - q * newt;
        t = newt;
        newt = tmp_t;
        const int64_t tmp_r = r - q * newr;
        r = newr;
        newr = tmp_r;
    }

    if (r != 1) {
        throw std::invalid_argument("ModInverseU32: no inverse exists");
    }

    if (t < 0) {
        t += static_cast<int64_t>(mod);
    }
    return static_cast<uint32_t>(t);
}

struct MPQSPolynomial {
    mpz_class a;
    mpz_class b;
    mpz_class c;
    std::vector<uint32_t> aPrimeIndices; // indices into factor base primes
};

struct FBPrime {
    uint32_t p;
    uint32_t r1;
    uint32_t r2;
};

struct Relation {
    mpz_class x; // x = sqrtN + offset
    std::bitset<kMaxCols> parity; // bit0 = sign, bits 1.. = primes
    std::vector<std::pair<uint32_t, uint16_t>> exps; // (primeIndex, exponent)
    std::vector<std::pair<uint64_t, uint16_t>> extraExps; // primes not in factor base (e.g. cancelled large primes)
};

static std::optional<std::pair<mpz_class, mpz_class>>
TryDependency(
    const mpz_class& N,
    const std::vector<uint32_t>& Primes,
    const std::vector<Relation>& Relations,
    const std::bitset<kMaxRelations>& comb,
    size_t relCount
) {
    mpz_class X = 1;
    std::vector<uint32_t> expSums(Primes.size(), 0);
    std::unordered_map<uint64_t, uint32_t> extraSums;

    for (size_t i = 0; i < relCount; ++i) {
        if (!comb.test(i)) {
            continue;
        }

        X = (X * Relations[i].x) % N;
        for (const auto& [primeIndex, exponent] : Relations[i].exps) {
            expSums[primeIndex] += exponent;
        }
        for (const auto& [primeValue, exponent] : Relations[i].extraExps) {
            extraSums[primeValue] += exponent;
        }
    }

    mpz_class Y = 1;
    for (size_t pi = 0; pi < Primes.size(); ++pi) {
        const uint32_t half = expSums[pi] / 2;
        if (half == 0) continue;

        mpz_class base = Primes[pi];
        mpz_class pow;
        mpz_powm_ui(pow.get_mpz_t(), base.get_mpz_t(), half, N.get_mpz_t());
        Y = (Y * pow) % N;
    }

    for (const auto& [primeValue, sumExp] : extraSums) {
        const uint32_t half = sumExp / 2;
        if (half == 0) {
            continue;
        }
        mpz_class base = primeValue;
        mpz_class pow;
        mpz_powm_ui(pow.get_mpz_t(), base.get_mpz_t(), half, N.get_mpz_t());
        Y = (Y * pow) % N;
    }

    mpz_class diff = X >= Y ? (X - Y) : (Y - X);
    mpz_class g1 = primetools::Gcd(diff, N);
    if (g1 != 1 && g1 != N) {
        return std::make_pair(g1, N / g1);
    }

    mpz_class sum = X + Y;
    mpz_class g2 = primetools::Gcd(sum, N);
    if (g2 != 1 && g2 != N) {
        return std::make_pair(g2, N / g2);
    }

    return std::nullopt;
}

static std::optional<std::pair<mpz_class, mpz_class>>
SolveRelations(
    const mpz_class& N,
    const std::vector<uint32_t>& Primes,
    const std::vector<Relation>& Relations,
    const size_t ColBits,
    LogCallback LogFn
) {
    const size_t relCount = Relations.size();
    if (ColBits > kMaxCols || relCount > kMaxRelations) {
        return std::nullopt;
    }

    const bool debug = (std::getenv("PRIMETOOLS_QS_DEBUG") != nullptr);

    struct Row {
        std::bitset<kMaxCols> cols;
        std::bitset<kMaxRelations> comb;
    };

    std::vector<Row> rows;
    rows.reserve(relCount);
    for (size_t i = 0; i < relCount; ++i) {
        Row r;
        r.cols = Relations[i].parity;
        r.comb.reset();
        r.comb.set(i);
        rows.push_back(std::move(r));
    }

    size_t pivotRow = 0;
    for (size_t col = 0; col < ColBits && pivotRow < rows.size(); ++col) {
        size_t best = pivotRow;
        while (best < rows.size() && !rows[best].cols.test(col)) {
            ++best;
        }
        if (best == rows.size()) {
            continue;
        }
        if (best != pivotRow) {
            std::swap(rows[best], rows[pivotRow]);
        }

        for (size_t r = 0; r < rows.size(); ++r) {
            if (r == pivotRow) continue;
            if (rows[r].cols.test(col)) {
                rows[r].cols ^= rows[pivotRow].cols;
                rows[r].comb ^= rows[pivotRow].comb;
            }
        }
        ++pivotRow;
    }

    for (const auto& row : rows) {
        if (row.cols.none()) {
            if (debug) {
                // We could have many dependencies; only log the first few attempts.
                static size_t depLogged = 0;
                if (depLogged < 3) {
                    LogFn("[MPQS] dependency candidate");
                    ++depLogged;
                }
            }
            auto f = TryDependency(N, Primes, Relations, row.comb, relCount);
            if (f) return f;
        }
    }

    if (debug) {
        size_t depCount = 0;
        for (const auto& row : rows) {
            if (row.cols.none()) {
                ++depCount;
            }
        }
        LogFn(std::string("[MPQS] dependencies=") + std::to_string(depCount));
    }

    return std::nullopt;
}

static void
MergeExtraExps(
    std::vector<std::pair<uint64_t, uint16_t>>& out,
    const std::vector<std::pair<uint64_t, uint16_t>>& a,
    const std::vector<std::pair<uint64_t, uint16_t>>& b
);

static Relation
CombineRelations(
    const mpz_class& N,
    const Relation& A,
    const Relation& B,
    const uint64_t CancelledLargePrime
)
{
    Relation out;
    out.x = (A.x * B.x) % N;
    out.parity = A.parity;
    out.parity ^= B.parity;

    // Merge exponent lists by prime index (both are naturally sorted).
    out.exps.reserve(A.exps.size() + B.exps.size());
    size_t ia = 0;
    size_t ib = 0;
    while (ia < A.exps.size() || ib < B.exps.size()) {
        if (ib >= B.exps.size() || (ia < A.exps.size() && A.exps[ia].first < B.exps[ib].first)) {
            out.exps.push_back(A.exps[ia++]);
        } else if (ia >= A.exps.size() || B.exps[ib].first < A.exps[ia].first) {
            out.exps.push_back(B.exps[ib++]);
        } else {
            // Same prime index.
            const uint32_t idx = A.exps[ia].first;
            const uint32_t sum = static_cast<uint32_t>(A.exps[ia].second) + static_cast<uint32_t>(B.exps[ib].second);
            out.exps.emplace_back(idx, static_cast<uint16_t>(std::min<uint32_t>(sum, std::numeric_limits<uint16_t>::max())));
            ++ia;
            ++ib;
        }
    }

    // Merge any extra primes not represented in the factor base.
    MergeExtraExps(out.extraExps, A.extraExps, B.extraExps);

    // If we combined two partial relations with the same large prime, that prime
    // contributes a square (p^2) to the product, so it must be included in Y.
    if (CancelledLargePrime != 0) {
        out.extraExps.emplace_back(CancelledLargePrime, 2);
        if (out.extraExps.size() > 1) {
            std::sort(out.extraExps.begin(), out.extraExps.end(),
                      [](const auto& x, const auto& y) { return x.first < y.first; });
            size_t w = 0;
            for (size_t i = 0; i < out.extraExps.size(); ++i) {
                if (w == 0 || out.extraExps[i].first != out.extraExps[w - 1].first) {
                    out.extraExps[w++] = out.extraExps[i];
                } else {
                    const uint32_t sum = static_cast<uint32_t>(out.extraExps[w - 1].second) +
                                         static_cast<uint32_t>(out.extraExps[i].second);
                    out.extraExps[w - 1].second = static_cast<uint16_t>(
                        std::min<uint32_t>(sum, std::numeric_limits<uint16_t>::max()));
                }
            }
            out.extraExps.resize(w);
        }
    }

    return out;
}

static void
MergeExtraExps(
    std::vector<std::pair<uint64_t, uint16_t>>& out,
    const std::vector<std::pair<uint64_t, uint16_t>>& a,
    const std::vector<std::pair<uint64_t, uint16_t>>& b
)
{
    out.clear();
    out.reserve(a.size() + b.size() + 1);
    out.insert(out.end(), a.begin(), a.end());
    out.insert(out.end(), b.begin(), b.end());
    if (out.size() > 1) {
        std::sort(out.begin(), out.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
        size_t w = 0;
        for (size_t i = 0; i < out.size(); ++i) {
            if (w == 0 || out[i].first != out[w - 1].first) {
                out[w++] = out[i];
            } else {
                const uint32_t sum = static_cast<uint32_t>(out[w - 1].second) + static_cast<uint32_t>(out[i].second);
                out[w - 1].second = static_cast<uint16_t>(std::min<uint32_t>(sum, std::numeric_limits<uint16_t>::max()));
            }
        }
        out.resize(w);
    }
}

// Attempt to factor Q(x) = x^2 - N over the factor base.
// Returns:
//  - true + LargePrime==0: fully smooth relation
//  - true + LargePrime>0: partial relation with one leftover prime LargePrime
//  - false: not usable
static bool
BuildRelation(
    const mpz_class& N,
    const mpz_class& X,
    const mpz_class& Value,
    const std::vector<uint32_t>& ExtraPrimeIndices,
    const std::vector<uint32_t>& Primes,
    const uint64_t LargePrimeBound,
    Relation& Out,
    uint64_t& LargePrime
) {
    mpz_class Q = Value;
    const bool negative = (Q < 0);
    if (negative) {
        Q = -Q;
    }

    if (Primes.size() + 1 > kMaxCols) {
        return false;
    }

    std::bitset<kMaxCols> parity;
    parity.reset();
    if (negative) {
        parity.set(0);
    }

    std::vector<std::pair<uint32_t, uint16_t>> exps;
    exps.reserve(32);

    for (uint32_t idx : ExtraPrimeIndices) {
        if (idx >= Primes.size()) {
            return false;
        }
        parity.flip(static_cast<size_t>(idx) + 1);
        exps.emplace_back(idx, 1);
    }

    for (uint32_t idx = 0; idx < Primes.size(); ++idx) {
        const uint32_t p = Primes[idx];
        uint16_t e = 0;
        while (mpz_divisible_ui_p(Q.get_mpz_t(), p)) {
            mpz_divexact_ui(Q.get_mpz_t(), Q.get_mpz_t(), p);
            ++e;
        }
        if (e) {
            if (e & 1u) {
                parity.flip(static_cast<size_t>(idx) + 1);
            }
            exps.emplace_back(idx, e);
        }
        if (Q == 1) {
            break;
        }
    }

    if (exps.size() > 1) {
        std::sort(exps.begin(), exps.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        size_t w = 0;
        for (size_t i = 0; i < exps.size(); ++i) {
            if (w == 0 || exps[i].first != exps[w - 1].first) {
                exps[w++] = exps[i];
            } else {
                const uint32_t sum = static_cast<uint32_t>(exps[w - 1].second) + static_cast<uint32_t>(exps[i].second);
                exps[w - 1].second = static_cast<uint16_t>(std::min<uint32_t>(sum, std::numeric_limits<uint16_t>::max()));
            }
        }
        exps.resize(w);
    }

    if (Q != 1) {
        // Single large-prime variant: accept if remaining cofactor is one prime <= bound.
        if (Q > 1 && mpz_fits_ulong_p(Q.get_mpz_t())) {
            const uint64_t lp = mpz_get_ui(Q.get_mpz_t());
            if (lp <= LargePrimeBound && primetools::IsPrime(lp)) {
                LargePrime = lp;
                Out.x = X;
                Out.parity = std::move(parity);
                Out.exps = std::move(exps);
                return true;
            }
        }
        return false;
    }

    Out.x = X;
    Out.parity = std::move(parity);
    Out.exps = std::move(exps);
    LargePrime = 0;
    return true;
}

static std::optional<MPQSPolynomial>
MakeMPQSPolynomial(
    const mpz_class& N,
    const mpz_class& sqrtN,
    const std::vector<uint32_t>& primes,
    const std::vector<FBPrime>& fb,
    const size_t polyIndex,
    const size_t aPrimeCount,
    std::mt19937_64& rng
)
{
    if (aPrimeCount == 0) {
        return std::nullopt;
    }
    if (primes.size() != fb.size()) {
        return std::nullopt;
    }

    // Choose primes for a so that a is near the MPQS target.
    // A common heuristic is aTarget ≈ 2*sqrt(N) / M, where M is the sieve interval.
    // We can't see M directly here, but we can infer a reasonable prime scale from sqrt(N)
    // and the caller-chosen sieve interval.
    //
    // In practice, we approximate the desired prime size as pTarget ≈ aTarget^(1/k).
    const double lnS = primetools::ln_mpz(sqrtN);
    if (!(lnS > 0.0)) {
        return std::nullopt;
    }

    // Caller uses M = sieveInterval and xOffset=M/2.
    // Empirically this implementation uses M in the 2^13..2^16 range.
    // Use the largest (2^16) as a conservative default for the target.
    const double lnM = std::log(static_cast<double>(1u << 16));
    const double lnATarget = std::max<double>(0.0, std::log(2.0) + lnS - lnM);
    const double pTargetD = std::exp(lnATarget / static_cast<double>(aPrimeCount));

    const uint32_t pMax = primes.back();
    const uint32_t pTarget = static_cast<uint32_t>(std::clamp<double>(pTargetD, 3.0, static_cast<double>(pMax)));
    uint32_t lowVal = std::max<uint32_t>(3u, pTarget / 2u);
    uint32_t highVal = std::min<uint32_t>(pMax, pTarget * 2u);
    if (highVal <= lowVal) {
        lowVal = 3u;
        highVal = pMax;
    }

    auto loIt = std::lower_bound(primes.begin(), primes.end(), lowVal);
    auto hiIt = std::upper_bound(primes.begin(), primes.end(), highVal);
    if (loIt == primes.end() || loIt >= hiIt) {
        loIt = primes.begin() + 1; // skip 2
        hiIt = primes.end();
    }

    const size_t loIdx = static_cast<size_t>(std::distance(primes.begin(), loIt));
    const size_t hiIdx = static_cast<size_t>(std::distance(primes.begin(), hiIt));
    if (hiIdx <= loIdx + 1) {
        return std::nullopt;
    }

    std::uniform_int_distribution<size_t> dist(loIdx, hiIdx - 1);

    std::unordered_set<uint32_t> chosen;
    chosen.reserve(aPrimeCount * 2);

    std::vector<uint32_t> aPrimeIndices;
    aPrimeIndices.reserve(aPrimeCount);

    for (size_t tries = 0; tries < 256 && aPrimeIndices.size() < aPrimeCount; ++tries) {
        const uint32_t idx = static_cast<uint32_t>(dist(rng));
        if (idx == 0) {
            continue;
        }
        const uint32_t p = primes[idx];
        if (p <= 2) {
            continue;
        }
        if (chosen.insert(p).second) {
            aPrimeIndices.push_back(idx);
        }
    }
    if (aPrimeIndices.size() != aPrimeCount) {
        return std::nullopt;
    }

    // CRT to find b0 such that b0 ≡ ±sqrt(N) (mod p_i) for each p_i|a.
    mpz_class b0 = 0;
    mpz_class m = 1;

    for (size_t j = 0; j < aPrimeIndices.size(); ++j) {
        const uint32_t idx = aPrimeIndices[j];
        const uint32_t p = primes[idx];
        const uint32_t r1 = fb[idx].r1;
        const uint32_t r2 = fb[idx].r2;
        const uint32_t residue = (((polyIndex >> j) & 1u) && r2 != r1) ? r2 : r1;

        const uint32_t m_mod_p = static_cast<uint32_t>(primetools::modulo(m, p));
        if (m_mod_p == 0) {
            return std::nullopt;
        }
        uint32_t inv;
        try {
            inv = ModInverseU32(m_mod_p, p);
        } catch (...) {
            return std::nullopt;
        }
        const uint32_t b0_mod_p = static_cast<uint32_t>(primetools::modulo(b0, p));
        const uint32_t diff = (residue >= b0_mod_p) ? (residue - b0_mod_p) : (residue + p - b0_mod_p);
        const uint32_t t = static_cast<uint32_t>((static_cast<uint64_t>(diff) * inv) % p);

        b0 += m * static_cast<unsigned long>(t);
        m *= static_cast<unsigned long>(p);
        b0 %= m;
    }

    MPQSPolynomial poly;
    poly.a = m;
    poly.aPrimeIndices = std::move(aPrimeIndices);

    // Choose b close to sqrt(N) to keep coefficients well-behaved.
    mpz_class b = b0;
    if (poly.a != 0) {
        mpz_class k = (sqrtN - b0) / poly.a;
        b = b0 + k * poly.a;

        // Round to nearest multiple.
        mpz_class d = sqrtN - b;
        mpz_class half = poly.a / 2;
        if (d > half) {
            b += poly.a;
        } else if (-d > half) {
            b -= poly.a;
        }
    }
    poly.b = b;

    // c = (b^2 - N) / a
    mpz_class bb = poly.b * poly.b - N;
    if (!mpz_divisible_p(bb.get_mpz_t(), poly.a.get_mpz_t())) {
        return std::nullopt;
    }
    mpz_divexact(poly.c.get_mpz_t(), bb.get_mpz_t(), poly.a.get_mpz_t());

    return poly;
}

} // namespace

const uint32_t
QuadraticSieveRecommendedFactorBaseBound(
    const mpz_class& N
) {
    if (N < 2) {
        return 0;
    }

    const size_t bits = mpz_sizeinbase(N.get_mpz_t(), 2);

    // Practical piecewise heuristic for this implementation.
    // (The asymptotic QS estimate tends to undershoot for small N and explode for huge N.)
    if (bits <= 32) return 200;
    if (bits <= 48) return 500;
    if (bits <= 64) return 2000;
    if (bits <= 80) return 5000;
    if (bits <= 96) return 12000;
    if (bits <= 112) return 25000;
    if (bits <= 128) return 45000;

    const double lnN = primetools::ln_mpz(N);
    if (lnN <= 1.0) {
        return 5000;
    }
    const double lnlnN = std::log(lnN);

    // Standard QS-style form: exp(0.5 * sqrt(lnN * lnlnN))
    // Clamp to keep runtime/memory bounded for this basic implementation.
    const double t = 0.5 * std::sqrt(std::max(0.0, lnN * lnlnN));
    const double raw = std::exp(t);
    const uint64_t rounded = static_cast<uint64_t>(raw + 0.5);
    return static_cast<uint32_t>(std::clamp<uint64_t>(rounded, 50000u, 200000u));
}

const std::optional<std::pair<mpz_class, mpz_class>>
QuadraticSieveFactor(
    const mpz_class& N,
    LogCallback LogFn
) {
    if (N < 2) {
        return std::nullopt;
    }
    if ((N & 1) == 0) {
        return std::make_pair(mpz_class(2), N / 2);
    }
    if (primetools::IsPrime(N)) {
        return std::nullopt;
    }

    const mpz_class sqrtN = primetools::Sqrt(N);
    if (sqrtN * sqrtN == N) {
        return std::make_pair(sqrtN, sqrtN);
    }

    const size_t nBits = mpz_sizeinbase(N.get_mpz_t(), 2);
    const bool debug = (std::getenv("PRIMETOOLS_QS_DEBUG") != nullptr);

    // Choose a factor base bound. For small N we bias upward for reliability.
    uint32_t B = QuadraticSieveRecommendedFactorBaseBound(N);
    if (nBits <= 40) {
        B = std::max<uint32_t>(B, 1000u);
    } else if (nBits <= 64) {
        B = std::max<uint32_t>(B, 2000u);
    }

    // Choose sieving parameters based on size.
    size_t sieveInterval = 1u << 15;
    size_t extraRelations = 24;
    size_t maxBlocks = 96;

    if (nBits <= 48) {
        sieveInterval = 1u << 13;
        extraRelations = 24;
        maxBlocks = 96;
    } else if (nBits <= 80) {
        sieveInterval = 1u << 14;
        extraRelations = 28;
        maxBlocks = 128;
    } else if (nBits <= 112) {
        sieveInterval = 1u << 15;
        extraRelations = 32;
        maxBlocks = 160;
    } else {
        sieveInterval = 1u << 16;
        extraRelations = 96;
        maxBlocks = 384;
    }

    // Build factor base.
    std::vector<FBPrime> fb;
    std::vector<uint32_t> primes;
    fb.reserve(std::max<uint32_t>(B, 10u) / 10);
    primes.reserve(std::max<uint32_t>(B, 10u) / 10);

    // Always include p=2 with roots for odd N: x^2 = N (mod 2) => x = 1.
    fb.push_back(FBPrime{2u, 1u, 1u});
    primes.push_back(2u);

    PrimeGenerator<uint64_t> gen;
    for (;;) {
        const uint64_t p64 = gen.Next();
        if (p64 <= 2) {
            continue;
        }
        if (p64 > B) {
            break;
        }
        const uint32_t p = static_cast<uint32_t>(p64);

        const uint64_t nmodp = primetools::modulo(N, p);
        if (nmodp == 0) {
            return std::make_pair(mpz_class(p), N / p);
        }
        if (primetools::legendre_symbol(nmodp, p) != 1) {
            continue;
        }

        auto r = primetools::TonelliShanks<uint64_t>(nmodp, p);
        if (!r) {
            continue;
        }
        const uint32_t r1 = static_cast<uint32_t>(*r);
        const uint32_t r2 = (r1 == 0) ? 0u : (p - r1);

        fb.push_back(FBPrime{p, r1, r2});
        primes.push_back(p);
    }

    const size_t colBits = primes.size() + 1; // sign + primes
    const size_t targetRelations = colBits + extraRelations;

    const size_t candidatesPerBlock = [&]() -> size_t {
        if (nBits <= 64) return 2048;
        if (nBits <= 96) return 4096;
        if (nBits <= 128) return 8192;
        if (nBits <= 160) return 32768;
        return 49152;
    }();

    const size_t relationSlack = [&]() -> size_t {
        // Extra relations beyond the column count improves the odds of finding
        // a usable dependency (especially without large-prime variants).
        if (nBits <= 96) return 256;
        if (nBits <= 128) return 1024;
        if (nBits <= 160) return 4096;
        return 8192;
    }();

    const uint64_t largePrimeBound = [&]() -> uint64_t {
        // Single large-prime variant bound.
        // For this implementation, we intentionally cap it to encourage collisions
        // (pairing partial relations) rather than collecting mostly-unique large primes.
        const unsigned __int128 bb = static_cast<unsigned __int128>(B) * static_cast<unsigned __int128>(B);
        const unsigned __int128 bound = bb * 4u;
        const uint64_t cap = 10'000'000ULL;
        const unsigned __int128 maxu64 = static_cast<unsigned __int128>(std::numeric_limits<uint64_t>::max());
        const uint64_t unclamped = static_cast<uint64_t>(std::min(bound, maxu64));
        return std::min<uint64_t>(unclamped, cap);
    }();

    std::vector<Relation> relations;
    relations.reserve(targetRelations + 32);

    // Single-large-prime partial relation store: when we see the same large prime twice,
    // multiply the relations to cancel it and create a full relation.
    std::unordered_map<uint64_t, Relation> partialByLargePrime;
    std::mutex partialMutex;

    size_t nextSolveAt = targetRelations;

    // MPQS: iterate over multiple polynomials
    const size_t M = sieveInterval;
    const size_t xOffset = M / 2;
    std::mt19937_64 rng(static_cast<uint64_t>(mpz_get_ui(N.get_mpz_t())) ^ (static_cast<uint64_t>(nBits) << 32));
    const size_t aPrimeCount = [&]() -> size_t {
        if (nBits <= 96) return 3;
        if (nBits <= 160) return 4;
        return 5;
    }();

    for (size_t polyIdx = 0; polyIdx < maxBlocks; ++polyIdx) {
        auto polyOpt = MakeMPQSPolynomial(N, sqrtN, primes, fb, polyIdx, aPrimeCount, rng);
        if (!polyOpt) {
            continue;
        }
        const MPQSPolynomial& poly = *polyOpt;

        std::vector<double> score(M, 0.0);

        // Compute f(x) = a*x^2 + 2*b*x + c incrementally for log(|f|).
        const int64_t x0 = -static_cast<int64_t>(xOffset);
        mpz_class x = static_cast<long>(x0);
        mpz_class F = poly.a * x * x + 2 * poly.b * x + poly.c;
        mpz_class delta = poly.a * (2 * x + 1) + 2 * poly.b;
        const mpz_class deltaInc = 2 * poly.a;

        for (size_t i = 0; i < M; ++i) {
            mpz_class absF = (F < 0) ? -F : F;
            const size_t bits = mpz_sizeinbase(absF.get_mpz_t(), 2);
            score[i] = static_cast<double>(bits) * std::log(2.0);

            F += delta;
            delta += deltaInc;
        }

        // Sieve using modular roots for this polynomial.
        for (size_t pi = 0; pi < fb.size(); ++pi) {
            const uint32_t p = fb[pi].p;
            const double lp = std::log(static_cast<double>(p));

            const uint32_t a_mod = static_cast<uint32_t>(primetools::modulo(poly.a, p));
            const uint32_t b_mod = static_cast<uint32_t>(primetools::modulo(poly.b, p));
            const uint32_t xoff_mod = static_cast<uint32_t>(xOffset % p);

            auto sieve_root = [&](uint32_t xr) {
                const uint32_t offset = (xr + xoff_mod) % p;
                for (size_t j = offset; j < M; j += p) {
                    score[j] -= lp;
                }
            };

            if (a_mod != 0) {
                uint32_t inva;
                try {
                    inva = ModInverseU32(a_mod, p);
                } catch (...) {
                    continue;
                }

                const uint32_t r1 = fb[pi].r1;
                const uint32_t r2 = fb[pi].r2;

                const uint32_t d1 = (r1 >= b_mod) ? (r1 - b_mod) : (r1 + p - b_mod);
                const uint32_t x1 = static_cast<uint32_t>((static_cast<uint64_t>(d1) * inva) % p);
                sieve_root(x1);

                if (r2 != r1) {
                    const uint32_t d2 = (r2 >= b_mod) ? (r2 - b_mod) : (r2 + p - b_mod);
                    const uint32_t x2 = static_cast<uint32_t>((static_cast<uint64_t>(d2) * inva) % p);
                    if (x2 != x1) {
                        sieve_root(x2);
                    }
                }
            } else {
                // Special prime p | a: f(x) mod p is linear: 2*b*x + c == 0 (mod p)
                const uint32_t two_b = static_cast<uint32_t>((2ull * b_mod) % p);
                if (two_b == 0) {
                    continue;
                }
                uint32_t inv;
                try {
                    inv = ModInverseU32(two_b, p);
                } catch (...) {
                    continue;
                }
                const uint32_t c_mod = static_cast<uint32_t>(primetools::modulo(poly.c, p));
                const uint32_t rhs = (c_mod == 0) ? 0u : (p - c_mod);
                const uint32_t xr = static_cast<uint32_t>((static_cast<uint64_t>(rhs) * inv) % p);
                sieve_root(xr);
            }
        }

        // Candidate filtering: select the best-K candidates by sieve score.
        // This is more robust than a fixed threshold (we do not account for prime powers in the log sieve).
        const size_t K = std::min(candidatesPerBlock, M);
        std::vector<size_t> indices(M);
        std::iota(indices.begin(), indices.end(), 0);

        auto cmp = [&](size_t a, size_t b) { return score[a] < score[b]; };
        if (K < M) {
            std::nth_element(indices.begin(), indices.begin() + static_cast<ptrdiff_t>(K), indices.end(), cmp);
            indices.resize(K);
        }
        std::sort(indices.begin(), indices.end(), cmp);

        // Multi-threaded relation testing (dominant cost for large N).
        const size_t relationLimit = targetRelations + relationSlack;
        const size_t hwThreads = std::max<size_t>(1, std::thread::hardware_concurrency());
        const size_t threads = (indices.size() < 512) ? 1 : hwThreads;

        std::atomic<bool> stop{false};
        std::mutex relMutex;

        auto flush_local = [&](std::vector<Relation>& local) {
            if (local.empty()) {
                return;
            }
            std::scoped_lock lock(relMutex);
            for (auto& r : local) {
                if (relations.size() >= relationLimit) {
                    stop.store(true, std::memory_order_relaxed);
                    break;
                }
                relations.push_back(std::move(r));
            }
            local.clear();
            if (relations.size() >= relationLimit) {
                stop.store(true, std::memory_order_relaxed);
            }
        };

        auto worker = [&](size_t start, size_t end) {
            std::vector<Relation> local;
            local.reserve(16);

            for (size_t k = start; k < end; ++k) {
                if (stop.load(std::memory_order_relaxed)) {
                    break;
                }

                const size_t i = indices[k];
                Relation rel;
                const int64_t xi = static_cast<int64_t>(i) - static_cast<int64_t>(xOffset);
                const mpz_class xmp = static_cast<long>(xi);

                mpz_class X = poly.a * xmp + poly.b;
                X %= N;
                if (X < 0) {
                    X += N;
                }

                const mpz_class fval = poly.a * xmp * xmp + 2 * poly.b * xmp + poly.c;
                uint64_t largePrime = 0;
                if (!BuildRelation(N, X, fval, poly.aPrimeIndices, primes, largePrimeBound, rel, largePrime)) {
                    continue;
                }

                if (largePrime == 0) {
                    local.push_back(std::move(rel));
                } else {
                    std::optional<Relation> combined;
                    {
                        std::scoped_lock lock(partialMutex);
                        auto it = partialByLargePrime.find(largePrime);
                        if (it == partialByLargePrime.end()) {
                            partialByLargePrime.emplace(largePrime, std::move(rel));
                        } else {
                            Relation other = std::move(it->second);
                            partialByLargePrime.erase(it);
                            combined = CombineRelations(N, other, rel, largePrime);
                        }
                    }
                    if (combined) {
                        local.push_back(std::move(*combined));
                    }
                }

                if (local.size() >= 8) {
                    flush_local(local);
                }
            }

            flush_local(local);
        };

        if (threads == 1) {
            worker(0, indices.size());
        } else {
            std::vector<std::thread> pool;
            pool.reserve(threads);
            const size_t chunk = (indices.size() + threads - 1) / threads;

            for (size_t t = 0; t < threads; ++t) {
                const size_t start = t * chunk;
                const size_t end = std::min(indices.size(), start + chunk);
                if (start >= end) {
                    break;
                }
                pool.emplace_back(worker, start, end);
            }
            for (auto& th : pool) {
                th.join();
            }
        }

        if (relations.size() >= nextSolveAt) {
            auto factors = SolveRelations(N, primes, relations, colBits, LogFn);
            if (factors) {
                return factors;
            }
            // Try again after we have gathered a bit more.
            nextSolveAt = relations.size() + 256;
        }

        if (debug) {
            LogFn(std::string("[MPQS] poly=") + std::to_string(polyIdx) +
                  " a=" + poly.a.get_str() +
                  " relations=" + std::to_string(relations.size()) +
                  " partial=" + std::to_string(partialByLargePrime.size()));
        }
    }

    return std::nullopt;
}

} // namespace primetools
