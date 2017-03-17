
#include "MultiComposeMod.h"
#include "Util.h"
#include <thread>

using namespace std;

MultiComposeMod::MultiComposeMod() {
}

MultiComposeMod::~MultiComposeMod() {
}

void MultiComposeMod::compose(Vec<ZZ_pX>& result,
                              const Vec<ZZ_pX> &gVec,
                              const ZZ_pX &h,
                              const ZZ_pXModulus &F) {

    Mat<ZZ_p> hMat;
    ZZ_pXMultiplier hMultiplier;

    precompute(hMat, hMultiplier, gVec.length(), h, F);
    compose(result, gVec, hMat, hMultiplier, F);

    hMat.kill();
}

void MultiComposeMod::compose(ZZ_pX& result,
                              const ZZ_pX &g,
                              const ZZ_pX &h,
                              const ZZ_pXModulus &F) {
    CompMod(result, g, h, F);
//    Vec<ZZ_pX> tempVec;
//    tempVec.SetLength(1);
//    tempVec[0] = g;
//    compose(tempVec, tempVec, h, F);
//    result = tempVec[0];
//    tempVec.kill();
}

void MultiComposeMod::compose(ZZ_pX& result,
                              const ZZ_pX &g,
                              const Mat<ZZ_p> &hMat,
                              const ZZ_pXMultiplier &hMultiplier,
                              const ZZ_pXModulus &F) {
    Vec<ZZ_pX> tempVec;
    tempVec.SetLength(1);
    tempVec[0] = g;
    compose(tempVec, tempVec, hMat, hMultiplier, F);
    result = tempVec[0];
    tempVec.kill();
}


void MultiComposeMod::compose(Vec<ZZ_pX>& result,
                              const Vec<ZZ_pX> &gVec,
                              const Mat<ZZ_p> &hMat,
                              const ZZ_pXMultiplier &hMultiplier,
                              const ZZ_pXModulus &F) {

    long n = deg(F);
    long k = gVec.length();
    long t = SqrRoot(n * k);

    // number of subpolynomials in each polynomial
    long numSubPolys = (n + t - 1) / t;

    long gMatNumRows = k * numSubPolys;
    Mat<ZZ_p> gMat;
    gMat.SetDims(gMatNumRows, t);
    
    long rowIndex = 0;
    long colIndex = 0;
    
    // putting the g's in a matrix
    for (long i = 0; i < k; i++) {
        for (long j = 0; j < numSubPolys; j++) {
            for (long l = 0; l < t; l++) {
                gMat[rowIndex + j][l] = coeff(gVec[i], colIndex + l);
            }
            
            colIndex += t;
        }
        
        rowIndex += numSubPolys;
        colIndex = 0;
    }

    Mat<ZZ_p> resultMat;
    resultMat.SetDims(gMatNumRows, n);
    
    mul(resultMat, gMat, hMat);
    
    ZZ_pX tempG;
    clear(tempG);

    rowIndex = 0;
    for (long i = 0; i < k; i++) {

        // starting from the highest degree to use in Horner evaluation
        conv(result[i], resultMat[rowIndex + numSubPolys - 1]);
        for (long j = numSubPolys - 2; j >= 0; j--) {
            conv(tempG, resultMat[rowIndex + j]);
            MulMod(result[i], result[i], hMultiplier, F);
            add(result[i], result[i], tempG);
        }
        
        rowIndex += numSubPolys;
    }
    
    gMat.kill();
    tempG.kill();
}

void MultiComposeMod::precompute(Mat<ZZ_p> &hMat,
                                 ZZ_pXMultiplier &hMultiplier,
                                 long k,
                                 const ZZ_pX &h,
                                 const ZZ_pXModulus &F) {
    long n = deg(F);
    long t = SqrRoot(n * k);

    ZZ_pX tempH;
    set(tempH);

    hMat.SetDims(t, n);

    ZZ_pXMultiplier tempMultiplier;
    build(tempMultiplier, h, F);
    
    // putting the h^i in a matrix
    // after this loop we have tempH = h^t mod F
    for (long i = 0; i < t; i++) {
        VectorCopy(hMat[i], tempH, n);
        MulMod(tempH, tempH, tempMultiplier, F);
    }

    build(hMultiplier, tempH, F);
    tempH.kill();
}