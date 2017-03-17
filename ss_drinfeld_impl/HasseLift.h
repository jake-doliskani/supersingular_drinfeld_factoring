
#ifndef HASSELIFT_H
#define HASSELIFT_H

#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

using namespace NTL;

class HasseLift {
public:

    HasseLift(const ZZ_pX &g,
              const ZZ_pX &delta,
              const ZZ_pXModulus &F);

    void compute(ZZ_pX& result, long n);

    virtual ~HasseLift();

private:

    void computeFrobProduct(Mat<ZZ_pX> &result,
                            const Vec<Mat<ZZ_pX>> &B,
                            long l);

    void evaluateMatrices(Vec<Mat<ZZ_pX>> &result,
                          const Vec<ZZ_pX> &inverseFrobes,
                          long l);

    const double BETA = 0.8;
    ZZ_pX frobenius;
    ZZ_pX g;
    ZZ_pX delta;
    ZZ_pXModulus F;
};

#endif /* HASSELIFT_H */

