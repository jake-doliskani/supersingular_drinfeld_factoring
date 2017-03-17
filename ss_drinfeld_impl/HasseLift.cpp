#include "HasseLift.h"
#include <cmath>
#include "Util.h"
#include <iostream>
#include "BaseChange.h"
#include "MultiComposeMod.h"
#include "FrobComp.h"
#include "BalancedMul.h"


using namespace std;

HasseLift::HasseLift(const ZZ_pX &g,
                     const ZZ_pX &delta,
                     const ZZ_pXModulus &F) {
    this->g = g;
    this->delta = delta;
    this->F = F;
    PowerXMod(this->frobenius, ZZ_p::modulus(), F);
}

HasseLift::~HasseLift() {
    this->g.kill();
    this->delta.kill();
    this->F.f.kill();
    this->frobenius.kill();
}

void HasseLift::compute(ZZ_pX& result, long n) {
    long l = (long) pow(n, BETA);
    long m = n / l;

    Util util;
    long start = util.getTimeMillis();
    // compute x, tau^{-l}(x), tau^{-2l}(x), ..., tau^{-lm}(x)
    Vec<ZZ_pX> inverseFrobs;
    inverseFrobs.SetLength(m + 1);
    FrobComp frobComp;
    frobComp.computeInverseFrobes(inverseFrobs, m, l, frobenius, F);
    cout << "computeInverseFrobes: " << util.getTimeMillis() - start << endl;

    start = util.getTimeMillis();
    // compute B(inverseFrobs[0]), ..., B(inverseFrobs[m])
    Vec<Mat < ZZ_pX>> BEval;
    BEval.SetLength(inverseFrobs.length());

    for (long i = 0; i < BEval.length(); i++)
        BEval[i].SetDims(2, 2);

    evaluateMatrices(BEval, inverseFrobs, l);
    cout << "evaluate BEval: " << util.getTimeMillis() - start << endl;

    start = util.getTimeMillis();
    // compute the product tau^{lm}(BEval(m))...tau^{l}(BEval(1)) BEval(0)
    Mat<ZZ_pX> frobProduct;
    frobProduct.SetDims(2, 2);
    computeFrobProduct(frobProduct, BEval, l);
    cout << "computeFrobProduct2: " << util.getTimeMillis() - start << endl;

    MulMod(result, frobProduct[1][1], g, F);
    add(result, result, frobProduct[1][0]);

    inverseFrobs.kill();
    BEval.kill();
    frobProduct.kill();
}

void HasseLift::computeFrobProduct(Mat<ZZ_pX> &result,
                                   const Vec<Mat<ZZ_pX>> &B,
                                   long l) {
    long m = B.length();

    ZZ_pX frob;
    ZZ_pX xq;
    SetX(xq);
    FrobComp frobComp;
    frobComp.computeFrobPower(xq, xq, l, frobenius, F);
    frob = xq;

    Mat<ZZ_pX> tempMat;
    tempMat.SetDims(2, 2);

    result = B[0];
    MultiComposeMod multiComposeMod;
    Mat<ZZ_p> xqMat;
    ZZ_pXMultiplier xqMultiplier;
    multiComposeMod.precompute(xqMat, xqMultiplier, 1, xq, F);
    Vec<ZZ_pX> tempVec;
    tempVec.SetLength(4);

    BalancedMul balancedMul;
    
    for (long i = 1; i < m; i++) {

        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++) {
                tempVec[2 * j + k] = B[i][j][k];
            }

        multiComposeMod.compose(tempVec, tempVec, frob, F);

        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++) {
                tempMat[j][k] = tempVec[2 * j + k];
            }

        balancedMul.mul_ZZ_pXMat(result, tempMat, result, F);
        multiComposeMod.compose(frob, frob, xqMat, xqMultiplier, F);
    }

    frob.kill();
    tempMat.kill();
    xqMat.kill();
    tempVec.kill();
}

void HasseLift::evaluateMatrices(Vec<Mat<ZZ_pX>> &result,
                                 const Vec<ZZ_pX> &inverseFrobes,
                                 long l) {

    long m = inverseFrobes.length();
    Vec<ZZ_pX> frobCompos;
    frobCompos.SetLength(3);
    frobCompos[0] = delta;
    frobCompos[1] = g;
    Vec<ZZ_pX> tempVec;
    tempVec.SetLength(2);

    MultiComposeMod multiComposeMod;
    multiComposeMod.compose(frobCompos[1], frobCompos[1], frobenius, F);
    MulMod(frobCompos[2], frobenius, frobCompos[0], F);
    NTL::negate(frobCompos[2], frobCompos[2]);

    // initialize the matrices to identity
    for (long i = 0; i < m; i++) {
        set(result[i][0][0]);
        clear(result[i][0][1]);
        clear(result[i][1][0]);
        set(result[i][1][1]);
    }


    Mat<ZZ_p> frobMat;
    ZZ_pXMultiplier frobMultiplier;
    multiComposeMod.precompute(frobMat, frobMultiplier, 3, frobenius, F);

    Vec<ZZ_pXMultiplier> iFrobMultipliers;
    iFrobMultipliers.SetLength(inverseFrobes.length());
    for (long i = 0; i < inverseFrobes.length(); i++)
        build(iFrobMultipliers[i], inverseFrobes[i], F);

    BalancedMul balancedMul;
    for (long i = 0; i < l; i++) {

        tempVec[1] = frobCompos[1];
        for (long j = 0; j < m; j++) {
            MulMod(tempVec[0], frobCompos[0], iFrobMultipliers[j], F);
            add(tempVec[0], tempVec[0], frobCompos[2]);
            balancedMul.mul_ZZ_pXMatSpec(result[j], tempVec, result[j], F);
        }

        if (i < l - 1) {
            multiComposeMod.compose(frobCompos, frobCompos, frobMat,
                                    frobMultiplier, F);
        }
    }

    tempVec.kill();
    frobCompos.kill();
    frobMat.kill();
    frobMultiplier.b.kill();
    iFrobMultipliers.kill();
}




