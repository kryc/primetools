#include "fermat.hpp"

#include <algorithm>
#include <array>
#include <cinttypes>
#include <iostream>
#include <map>

#include <gmpxx.h>

#include "factors.hpp"
#include "util.hpp"

namespace primetools {

const FermatAlgorithm
GetFermatAlgorithmFromString(
    const std::string_view AlgorithmStr
)
{
    if (AlgorithmStr == "fermat" || AlgorithmStr == "standardfermat") {
        return AlgFermat;
    }
    else if (AlgorithmStr == "fermat2" || AlgorithmStr == "fermatalgorithm2") {
        return AlgFermat2;
    }
    else if (AlgorithmStr == "mffv4" || AlgorithmStr == "modifiedfermatv4") {
        return AlgModifiedFermatV4;
    }
    else if (AlgorithmStr == "fmmod20precomp") {
        return AlgFMMod20Precomp;
    }
    else {
        return AlgFermat;
    }
}

const std::string_view
FermatAlgorithmToString(
    const FermatAlgorithm Algorithm
)
{
    switch (Algorithm) {
        case AlgFermat:
            return "Standard Fermat";
        case AlgFermat2:
            return "Fermat Algorithm 2";
        case AlgModifiedFermatV4:
            return "Modified Fermat Factorisation V4";
        case AlgFMMod20Precomp:
            return "Fermat Mod 20 Precomputed";
        default:
            return "Unknown Algorithm";
    }
}

}