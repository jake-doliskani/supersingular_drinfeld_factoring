
#include <assert.h>

#include "BalancedMul.h"

BalancedMul::BalancedMul() {
}

BalancedMul::~BalancedMul() {
}

void BalancedMul::compute(Vec<ZZ_pEX> &vec) {
    long n = vec.length();

    while (n > 1) {
        for (long i = 0; i < n / 2; i++) {
            mul(vec[i], vec[2 * i], vec[2 * i + 1]);
        }

        if (n % 2 == 1) {
            mul(vec[n / 2 - 1], vec[n / 2 - 1], vec[n - 1]);
        }

        n /= 2;
    }
}

void BalancedMul::compute(Vec<Mat<ZZ_pEX>>&vec,
                          const ZZ_pEXModulus& F) {
    long n = vec.length();

    while (n > 1) {
        for (long i = 0; i < n / 2; i++) {
            mul_ZZ_pEXMat(vec[i], vec[2 * i], vec[2 * i + 1], F);
        }

        if (n % 2 == 1) {
            mul_ZZ_pEXMat(vec[n / 2 - 1], vec[n / 2 - 1], vec[n - 1], F);
        }

        n /= 2;
    }
}

void BalancedMul::mul_ZZ_pXMat(Mat<ZZ_pX>& result,
                               const Mat<ZZ_pX>& A,
                               const Mat<ZZ_pX>& B,
                               const ZZ_pXModulus &F) {
    Mat<ZZ_pX> tempMat;
    tempMat.SetDims(2, 2);
    ZZ_pX tempPoly;

    for (long i = 0; i < 2; i++) {
        for (long j = 0; j < 2; j++) {
            MulMod(tempMat[i][j], A[i][0], B[0][j], F);
            MulMod(tempPoly, A[i][1], B[1][j], F);
            add(tempMat[i][j], tempMat[i][j], tempPoly);
        }
    }

    result = tempMat;

    tempMat.kill();
    tempPoly.kill();
}

void BalancedMul::mul_ZZ_pEXMat(Mat<ZZ_pEX>& result,
                                const Mat<ZZ_pEX>& A,
                                const Mat<ZZ_pEX>& B,
                                const ZZ_pEXModulus &F) {
    Mat<ZZ_pEX> tempMat;
    tempMat.SetDims(2, 2);
    ZZ_pEX tempPoly;

    for (long i = 0; i < 2; i++) {
        for (long j = 0; j < 2; j++) {
            MulMod(tempMat[i][j], A[i][0], B[0][j], F);
            MulMod(tempPoly, A[i][1], B[1][j], F);
            add(tempMat[i][j], tempMat[i][j], tempPoly);
        }
    }

    result = tempMat;

    tempMat.kill();
    tempPoly.kill();
}

void BalancedMul::mul_ZZ_pXMatSpec(Mat<ZZ_pX> &result,
                                   const Vec<ZZ_pX> vec,
                                   const Mat<ZZ_pX> A,
                                   const ZZ_pXModulus &F) {
    Mat<ZZ_pX> tempMat;
    tempMat.SetDims(2, 2);
    tempMat[0][0] = A[1][0];
    tempMat[0][1] = A[1][1];

    ZZ_pX temp;
    MulMod(tempMat[1][0], vec[0], A[0][0], F);
    MulMod(temp, vec[1], A[1][0], F);
    add(tempMat[1][0], tempMat[1][0], temp);

    MulMod(tempMat[1][1], vec[0], A[0][1], F);
    MulMod(temp, vec[1], A[1][1], F);
    add(tempMat[1][1], tempMat[1][1], temp);

    result = tempMat;

    tempMat.kill();
    temp.kill();
}

void BalancedMul::mul_ZZ_pEXMatSpec(Mat<ZZ_pEX> &result,
                                   const Vec<ZZ_pEX> vec,
                                   const Mat<ZZ_pEX> A,
                                   const ZZ_pEXModulus &F) {
    Mat<ZZ_pEX> tempMat;
    tempMat.SetDims(2, 2);
    tempMat[0][0] = A[1][0];
    tempMat[0][1] = A[1][1];

    ZZ_pEX temp;
    MulMod(tempMat[1][0], vec[0], A[0][0], F);
    MulMod(temp, vec[1], A[1][0], F);
    add(tempMat[1][0], tempMat[1][0], temp);

    MulMod(tempMat[1][1], vec[0], A[0][1], F);
    MulMod(temp, vec[1], A[1][1], F);
    add(tempMat[1][1], tempMat[1][1], temp);

    result = tempMat;

    tempMat.kill();
    temp.kill();
}

