#ifndef PRIMEGENERATOR_HPP
#define PRIMEGENERATOR_HPP

#include <array>
#include <cstdint>

#include "util.hpp"
#include "wheel.hpp"

namespace primetools
{

template<typename T>
class PossiblePrimeGenerator
{
public:
    PossiblePrimeGenerator(
        const size_t Modulus = 510510,
        const T StartValue = 1
    ) 
    {
        Initialize(Modulus, StartValue);
    }

    const T& Next(
        void
    )
    {
        // Handle small primes below the wheel
        if (m_Count < m_SmallPrimesForWheel.size() && m_CurrentValue == 1)
        {
            m_Temp = m_SmallPrimesForWheel[m_Count++];
            return m_Temp;
        }
        m_Count++;
        const T gap = GetNextGap();
        m_CurrentValue += gap;
        return m_CurrentValue;
    }

    const T& Current(
        void
    ) const
    {
        return m_CurrentValue;
    }

private:
    void Initialize(
        const size_t Modulus,
        const T StartValue
    )
    {
        m_GapArray = GetWheelGapsForModulus(Modulus);
        m_SmallPrimesForWheel = GetPrimesForWheelModulus(Modulus);
        m_Gapword = m_GapArray[0];
        m_CurrentValue = StartValue;
        if ((m_CurrentValue % Modulus == 0))
        {
            m_CurrentValue += 1;
        }
        else if ((m_CurrentValue % Modulus) != 1)
        {
            // Walk forward to align with the wheel
            while ((m_CurrentValue % Modulus) != 1)
            {
                m_CurrentValue++;
            }
        }
    }

    const uint64_t
    GetNextGap(
        void
    )
    {
        uint64_t gap;
        do
        {
            gap = (m_Gapword & kGapMask);
            m_Gapword >>= kBitsPerWheelGap;
            m_GapIndex++;
            if (m_GapIndex >= kGapsPerWord)
            {
                m_GapIndex = 0;
                m_WheelIndex++;
                if (m_WheelIndex >= m_GapArray.size())
                {
                    m_WheelIndex = 0;
                }
                m_Gapword = m_GapArray[m_WheelIndex];
            }
        } while (gap == 0);
        return gap;
    }

    std::span<const uint64_t> m_GapArray;
    std::span<const uint64_t> m_SmallPrimesForWheel;
    uint64_t m_Gapword;
    size_t m_Count = 0;
    size_t m_WheelIndex = 0;
    size_t m_GapIndex = 0;
    T m_CurrentValue;
    T m_Temp;
};

template<typename T>
class PrimeGenerator : public PossiblePrimeGenerator<T>
{
public:
    const T& Next(
        void
    )
    {
        for(;;) {
            const T& candidate = PossiblePrimeGenerator<T>::Next();
            if (primetools::isprime(candidate)) {
                return candidate;
            }
        }
    }
    void Skip(
        const size_t Count
    )
    {
        for (size_t i = 0; i < Count; ++i)
        {
            Next();
        }
    }
};

} // namespace primetools

#endif // PRIMEGENERATOR_HPP