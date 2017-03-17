
#include "HasseLiftExt.h"
#include "Util.h"
#include "BaseChange.h"
#include "MultipointEval.h"
#include "FrobComp.h"
#include "MultiComposeMod.h"
#include "BalancedMul.h"
#include "HasseLift.h"
#include <NTL/matrix.h>
#include <cmath>
#include <NTL/ZZX.h>

HasseLiftExt::HasseLiftExt(const ZZ_pX &g,
                           const ZZ_pX &delta,
                           const ZZ_pXModulus &F) {
    this->g = g;
    this->delta = delta;
    this->F = F;
    PowerXMod(this->frobenius, ZZ_p::modulus(), F);
}

HasseLiftExt::~HasseLiftExt() {
    this->g.kill();
    this->delta.kill();
    this->F.f.kill();
    this->frobenius.kill();
    this->evalModulus.f.kill();
}

void HasseLiftExt::computeNaive(ZZ_pX& result, long n) {
    long l = (long) pow(n, BETA);
    long m = n / l;
    n = m * l + l + 1;

    ZZ_pX r0;
    ZZ_pX r1;
    ZZ_pX gTemp;
    ZZ_pX deltaTemp;
    ZZ_pX xq;
    ZZ_pX x;
    ZZ_pX temp1;
    ZZ_pX temp2;

    set(r0);
    r1 = g;
    PowerMod(gTemp, g, ZZ_p::modulus(), F);
    deltaTemp = delta;
    SetX(x);
    xq = this->frobenius;

    for (long i = 2; i <= n; i++) {
        sub(temp1, xq, x);
        MulMod(temp1, temp1, deltaTemp, F);
        MulMod(temp1, temp1, r0, F);
        MulMod(temp2, gTemp, r1, F);
        r0 = r1;
        sub(r1, temp2, temp1);

        PowerMod(gTemp, gTemp, ZZ_p::modulus(), F);
        PowerMod(deltaTemp, deltaTemp, ZZ_p::modulus(), F);
        PowerMod(xq, xq, ZZ_p::modulus(), F);
    }

    result = r1;

    r0.kill();
    r1.kill();
    gTemp.kill();
    deltaTemp.kill();
    xq.kill();
    x.kill();
    temp1.kill();
    temp2.kill();
}

void HasseLiftExt::compute(ZZ_pX& result, long n, int verbose) {

    long l = (long) pow(n, BETA);
    long m = n / l;

    ZZ_pE::init(F.f);

    Mat<ZZ_pEX> B;
    B.SetDims(2, 2);
    buildInitialMatrix(B);

    Util util;
    long start = util.getTimeMillis();
    // compute x, tau^{-l}(x), tau^{-2l}(x), ..., tau^{-lm}(x)
    Vec<ZZ_pX> inverseFrobs;
    inverseFrobs.SetLength(m + 1);
    FrobComp frobComp;
    frobComp.computeInverseFrobes(inverseFrobs, m, l, frobenius, F);
    buildEvalModulus(inverseFrobs);
    if (verbose == 1)
        cout << "computeInverseFrobes: " << util.getTimeMillis() - start << endl;

    start = util.getTimeMillis();
    // compute the product tau^{l - 1}(A)...tau(A)A
    computeFrobProduct(B, B, l);
    if (verbose == 1)
        cout << "computeFrobProduct: " << util.getTimeMillis() - start << endl;

    start = util.getTimeMillis();
    // compute B(inverseFrobs[0]), ..., B(inverseFrobs[m])
    Mat<Vec < ZZ_pX>> BEval;
    BEval.SetDims(2, 2);
    evaluate(BEval, B, inverseFrobs);
    if (verbose == 1)
        cout << "evaluate BEval: " << util.getTimeMillis() - start << endl;

    start = util.getTimeMillis();
    // compute the product tau^{lm}(BEval(i))...tau^{l}(BEval(i)) BEval(i)
    Mat<ZZ_pX> frobProduct;
    frobProduct.SetDims(2, 2);
    computeFrobProduct(frobProduct, BEval, l);
    if (verbose == 1)
        cout << "computeFrobProduct2: " << util.getTimeMillis() - start << endl;

    MulMod(result, frobProduct[1][1], g, F);
    add(result, result, frobProduct[1][0]);

    B.kill();
    inverseFrobs.kill();
    BEval.kill();
    frobProduct.kill();
}

void HasseLiftExt::buildEvalModulus(const Vec<ZZ_pX> &inverseFrobs) {
    long m = inverseFrobs.length();
    Vec<ZZ_pEX> tempVec;
    tempVec.SetLength(m);

    for (long i = 0; i < m; i++) {
        clear(tempVec[i]);
        SetCoeff(tempVec[i], 0, -to_ZZ_pE(inverseFrobs[i]));
        SetCoeff(tempVec[i], 1, 1);
    }

    BalancedMul balancedMul;
    balancedMul.compute(tempVec);
    evalModulus = tempVec[0];
    tempVec.kill();
}

void HasseLiftExt::evaluate(Mat<Vec<ZZ_pX> >& result,
                            const Mat<ZZ_pEX>& A,
                            const Vec<ZZ_pX>& points) {

    MultipointEval multiPointEval(points);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {

            result[i][j].SetLength(points.length());
            multiPointEval.eval(result[i][j], A[i][j]);
        }
    }
}

void HasseLiftExt::buildInitialMatrix(Mat<ZZ_pEX>& result) {
    clear(result[0][0]);
    set(result[0][1]);

    ZZ_pX temp;

    MulMod(temp, frobenius, delta, F);
    NTL::negate(temp, temp);
    clear(result[1][0]);
    SetCoeff(result[1][0], 0, to_ZZ_pE(temp));
    SetCoeff(result[1][0], 1, to_ZZ_pE(delta));

    MultiComposeMod().compose(temp, g, frobenius, F);
    clear(result[1][1]);
    SetCoeff(result[1][1], 0, to_ZZ_pE(temp));

    temp.kill();
}

//void HasseLiftExt::computeFrobProduct(Mat<ZZ_pEX>& result,
//                                      const Mat<ZZ_pEX>& A,
//                                      long l) {
//    Vec<ZZ_pX> frobCompos;
//    frobCompos.SetLength(3);
//    frobCompos[0] = coeff(A[1][0], 1)._ZZ_pE__rep;
//    frobCompos[1] = coeff(A[1][1], 0)._ZZ_pE__rep;
//    frobCompos[2] = coeff(A[1][0], 0)._ZZ_pE__rep;
//
//    Mat<ZZ_p> frobMat;
//    ZZ_pXMultiplier frobMultiplier;
//    MultiComposeMod multiComposeMod;
//    multiComposeMod.precompute(frobMat, frobMultiplier, 3, frobenius, F);
//
//    Vec<Mat<ZZ_pEX>> tempVec;
//    tempVec.SetLength(l);
//    
//    for (long i = 0; i < l; i++)
//        tempVec[i].SetDims(2, 2);
//    
//    tempVec[l - 1] = A;
//
//    for (long i = l - 2; i >= 0; i--) {
//        multiComposeMod.compose(frobCompos, frobCompos, frobMat,
//                                frobMultiplier, F);
//        
//        clear(tempVec[i][0][0]);
//        set(tempVec[i][0][1]);
//        clear(tempVec[i][1][0]);
//        clear(tempVec[i][1][1]);        
//        SetCoeff(tempVec[i][1][0], 1, to_ZZ_pE(frobCompos[0]));
//        SetCoeff(tempVec[i][1][1], 0, to_ZZ_pE(frobCompos[1]));
//        SetCoeff(tempVec[i][1][0], 0, to_ZZ_pE(frobCompos[2]));
//    }
//
//    BalancedMul balancedMul;
//    balancedMul.compute(tempVec, evalModulus);
//    result = tempVec[0];
//
//    frobCompos.kill();
//    frobMat.kill();
//    frobMultiplier.b.kill();
//    for (long i = 0; i < l; i++)
//        tempVec[i].kill();
//    tempVec.kill();
//}


//void HasseLiftExt::computeFrobProduct(Mat<ZZ_pEX>& result,
//                                      const Mat<ZZ_pEX>& A,
//                                      long l) {
//    Vec<ZZ_pX> frobCompos;
//    frobCompos.SetLength(3);
//    frobCompos[0] = coeff(A[1][0], 1)._ZZ_pE__rep;
//    frobCompos[1] = coeff(A[1][1], 0)._ZZ_pE__rep;
//    frobCompos[2] = coeff(A[1][0], 0)._ZZ_pE__rep;
//
//
//    Vec<Mat < ZZ_pEX>> tempVec;
//    tempVec.SetLength(l);
//
//    for (long i = 0; i < l; i++)
//        tempVec[i].SetDims(2, 2);
//
//    tempVec[l - 1] = A;
//
//    Util util;
//    long start = util.getTimeMillis();
//    for (long i = l - 2; i >= 0; i--) {
//        for (long j = 0; j < frobCompos.length(); j++)
//            PowerMod(frobCompos[j], frobCompos[j], ZZ_p::modulus(), F);
//
//        clear(tempVec[i][0][0]);
//        set(tempVec[i][0][1]);
//        clear(tempVec[i][1][0]);
//        clear(tempVec[i][1][1]);
//        SetCoeff(tempVec[i][1][0], 1, to_ZZ_pE(frobCompos[0]));
//        SetCoeff(tempVec[i][1][1], 0, to_ZZ_pE(frobCompos[1]));
//        SetCoeff(tempVec[i][1][0], 0, to_ZZ_pE(frobCompos[2]));
//    }
//    cout << "powering time: " << util.getTimeMillis() - start << endl;
//    
//    start = util.getTimeMillis();
//    BalancedMul balancedMul;
//    balancedMul.compute(tempVec, evalModulus);
//    result = tempVec[0];
//    cout << "balanced mul time: " << util.getTimeMillis() - start << endl;
//
//    frobCompos.kill();
//    for (long i = 0; i < l; i++)
//        tempVec[i].kill();
//    tempVec.kill();
//}

void HasseLiftExt::computeFrobProduct(Mat<ZZ_pEX>& result,
                                      const Mat<ZZ_pEX>& A,
                                      long l) {

    Vec<ZZ_pEX> tempVec;
    Mat<ZZ_pEX> tempMat;

    tempVec.SetLength(2);
    tempVec[0] = A[1][0];
    tempVec[1] = A[1][1];
    tempMat = A;

    BalancedMul balancedMul;
    for (long i = l - 2; i >= 0; i--) {

        SetCoeff(tempVec[0], 0, power(coeff(tempVec[0], 0), ZZ_p::modulus()));
        SetCoeff(tempVec[0], 1, power(coeff(tempVec[0], 1), ZZ_p::modulus()));
        SetCoeff(tempVec[1], 0, power(coeff(tempVec[1], 0), ZZ_p::modulus()));

        balancedMul.mul_ZZ_pEXMatSpec(tempMat, tempVec, tempMat, evalModulus);
    }

    result = tempMat;

    tempVec.kill();
    tempMat.kill();
}

void HasseLiftExt::computeFrobProduct(Mat<ZZ_pX> &result,
                                      Mat<Vec<ZZ_pX>> &B,
                                      long l) {
    long m = B[0][0].length();
    ZZ_pX xqInit;
    SetX(xqInit);
    FrobComp frobComp;
    frobComp.computeFrobPower(xqInit, xqInit, l, frobenius, F);

    // compute tau^2, tau^4, ..., tau^{2^k}
    // where k = floor(log m)
    Vec<ZZ_pX> frobPowers;
    frobPowers.SetLength(NumBits(m - 1));
    frobPowers[0] = xqInit;
    MultiComposeMod multiComposeMod;
    for (long i = 1; i < frobPowers.length(); i++)
        multiComposeMod.compose(frobPowers[i], frobPowers[i - 1],
                                frobPowers[i - 1], F);

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            actFrobPowers(B[i][j], frobPowers);

    Mat<ZZ_pX> tempMat;
    tempMat.SetDims(2, 2);

    result[0][0] = B[0][0][0];
    result[0][1] = B[0][1][0];
    result[1][0] = B[1][0][0];
    result[1][1] = B[1][1][0];

    BalancedMul balancedMul;
    for (long i = 1; i < m; i++) {
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++) {
                tempMat[j][k] = B[j][k][i];
            }

        balancedMul.mul_ZZ_pXMat(result, tempMat, result, F);
    }

    xqInit.kill();
    frobPowers.kill();
    tempMat.kill();
}

void HasseLiftExt::actFrobPowers(Vec<ZZ_pX> &vec,
                                 const Vec<ZZ_pX> &frobPowers) {
    long n = vec.length();
    long m = 1;
    long frobIndex = 0;

    Vec<ZZ_pX> tempVec;
    MultiComposeMod multiComposeMod;

    while (m < n) {
        for (long i = 0; i < n; i++) {
            if ((i & m) != 0)
                tempVec.append(vec[i]);
        }

        multiComposeMod.compose(tempVec, tempVec, frobPowers[frobIndex], F);

        long k = 0;
        for (long i = 0; i < n; i++) {
            if ((i & m) != 0) {
                vec[i] = tempVec[k];
                k++;
            }
        }

        tempVec.SetLength(0);
        m <<= 1;
        frobIndex++;
    }

    tempVec.kill();
}

