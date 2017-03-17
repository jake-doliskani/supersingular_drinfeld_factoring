
#ifndef SIMULCOMPOSEMOD_H
#define SIMULCOMPOSEMOD_H

#include <NTL/ZZ_pX.h>
#include <NTL/matrix.h>

using namespace NTL;

class MultiComposeMod {
public:

    MultiComposeMod();
    virtual ~MultiComposeMod();

    /**
     * Computes {@code gVec[i](h) mod F} for i = 0,...,m - 1 where
     * m = gVec.length().
     * @param result
     * @param gVec
     * @param h
     * @param F
     */
    void compose(Vec<ZZ_pX>& result,
                 const Vec<ZZ_pX> &gVec,
                 const ZZ_pX &h,
                 const ZZ_pXModulus &F);

    /**
     * Computes {@code gVec[i](h) mod F} for i = 0,...,m - 1 where
     * m = gVec.length(). The parameters {@code hMat, hMultiplier} are
     * precomputed using {@link #precompute}.
     * @param result
     * @param gVec
     * @param hMat
     * @param hMultiplier
     * @param F
     */
    void compose(Vec<ZZ_pX>& result,
                 const Vec<ZZ_pX> &gVec,
                 const Mat<ZZ_p> &hMat,
                 const ZZ_pXMultiplier &hMultiplier,
                 const ZZ_pXModulus &F);

    /**
     * Computes {@code g(h) mod F}.
     * @param result
     * @param g
     * @param h
     * @param F
     */
    void compose(ZZ_pX& result,
                 const ZZ_pX &g,
                 const ZZ_pX &h,
                 const ZZ_pXModulus &F);

    /**
     * Computes {@code g(h) mod F}. The parameters {@code hMat, hMultiplier} are
     * precomputed using {@link #precompute}.
     * @param result
     * @param g
     * @param hMat
     * @param hMultiplier
     * @param F
     */
    void compose(ZZ_pX& result,
                 const ZZ_pX &g,
                 const Mat<ZZ_p> &hMat,
                 const ZZ_pXMultiplier &hMultiplier,
                 const ZZ_pXModulus &F);

    /**
     * Pre-computes The parameters needed for modular composition. See, 
     * for example, {@link #compose}.
     * @param hMat
     * @param hMultiplier
     * @param k
     * @param h
     * @param F
     */
    void precompute(Mat<ZZ_p> &hMat,
                    ZZ_pXMultiplier &hMultiplier,
                    long k,
                    const ZZ_pX &h,
                    const ZZ_pXModulus &F);

private:

};

#endif /* SIMULCOMPOSEMOD_H */

