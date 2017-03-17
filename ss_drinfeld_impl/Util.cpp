
#include "Util.h"
#include <chrono>

using namespace std;

Util::Util() {
}

Util::~Util() {
}

long Util::getTimeMillis() {
    auto now = chrono::system_clock::now();
    auto now_ms = chrono::time_point_cast<chrono::milliseconds>(now);
    auto value = now_ms.time_since_epoch();
    long duration = value.count();
    return duration;
}

long Util::reverseBits(long a) {
    long result = 0;
    // reverse the bits of numPoints
    long temp = a;
    while (temp != 0) {
        result |= (temp & 0x1);
        result <<= 1;
        temp >>= 1;
    }

    result >>= 1;
    return result;
}

void Util::randomMonic(ZZ_pX& x, int degree) {
    random(x, degree);
    if (deg(x) < degree)
        SetCoeff(x, degree, 1);
}

