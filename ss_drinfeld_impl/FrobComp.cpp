
#include "FrobComp.h"
#include "BaseChange.h"
#include "MultiComposeMod.h"

FrobComp::FrobComp() {
}

FrobComp::~FrobComp() {
}

void FrobComp::computeInverseFrobes(Vec<ZZ_pX>& result,
                                    long m,
                                    long l,
                                    const ZZ_pX &initFrob,
                                    const ZZ_pXModulus &F) {
    ZZ_pX x;
    SetX(x);

    // compute x^{q^{-l}} using transposed modular composition
    ZZ_pX inverseFrob;
    computeFrobPower(inverseFrob, x, l, initFrob, F);
    BaseChange baseChange;
    baseChange.changeBasis(inverseFrob, inverseFrob, x, F);

    // compute the powers x, x^{q^{-l}}, x^{q^{-2l}}, ..., x^{q^{-ml}}
    computeFrobePowers(result, m, inverseFrob, F);

    // move the indices to make room for x
    for (long i = m; i > 0; i--) {
        result[i] = result[i - 1];
    }

    result[0] = x;

    x.kill();
    inverseFrob.kill();
}

void FrobComp::computeFrobePowers(Vec<ZZ_pX>& result,
                                  long m,
                                  const ZZ_pX &xqInit,
                                  const ZZ_pXModulus &F) {
    if (m == 1) {
        result[0] = xqInit;
        return;
    }

    MultiComposeMod multiComposeMod;
    
    if (m % 2 == 0) {

        Vec<ZZ_pX> temp;
        temp.SetLength(m / 2);
        computeFrobePowers(temp, m / 2, xqInit, F);

        for (long i = 0; i < m / 2; i++)
            result[i] = temp[i];

        multiComposeMod.compose(temp, temp, temp[m / 2 - 1], F);

        for (long i = m / 2; i < m; i++)
            result[i] = temp[i - m / 2];

        temp.kill();

    } else {

        computeFrobePowers(result, m - 1, xqInit, F);
        ZZ_pX temp;
        multiComposeMod.compose(temp, result[m - 2], xqInit, F);
        result[m - 1] = temp;
        temp.kill();
    }
}

void FrobComp::computeFrobPower(ZZ_pX& result,
                                const ZZ_pX &f,
                                long n,
                                const ZZ_pX &initFrob,
                                const ZZ_pXModulus &F) {
    if (n == 0) {
        result = f;
        return;
    }

    ZZ_pX xq;
    xq = initFrob;

    ZZ_pX temp;
    SetX(temp);

    MultiComposeMod multiComposeMod;
    
    long m = n;
    while (m != 0) {
        if ((m & 0x1) == 1)
            multiComposeMod.compose(temp, temp, xq, F);

        multiComposeMod.compose(xq, xq, xq, F);
        m >>= 1;
    }

    multiComposeMod.compose(result, f, temp, F);

    xq.kill();
    temp.kill();
}
