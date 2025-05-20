#ifndef EULER_HPP
#define EULER_HPP

template <typename T>
const T
EulerTotient(
    const T N
)
{
    T result = N;
    for (T p = 2; p * p <= N; ++p) {
        if (N % p == 0) {
            while (N % p == 0) {
                N /= p;
            }
            result -= result / p;
        }
    }
    if (N > 1) {
        result -= result / N;
    }
    return result;
}

#endif // EULER_HPP