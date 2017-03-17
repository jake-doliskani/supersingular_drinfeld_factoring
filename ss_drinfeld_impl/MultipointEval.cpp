
#include "MultipointEval.h"
#include "Util.h"

void MultipointEval::init(const Vec<ZZ_pE>& points) {
    this->numPoints = points.length();
    initSubProductTree(subProductTree);
    buildSubProductTree(points);
}

MultipointEval::MultipointEval(const Vec<ZZ_pE>& points) {
    init(points);
}

MultipointEval::MultipointEval(const Vec<ZZ_pX>& points) {
    long numPoints = points.length();
    
    Vec<ZZ_pE> tempPoints;
    tempPoints.SetLength(numPoints);

    for (long i = 0; i < numPoints; i++) {
        conv(tempPoints[i], points[i]);
    }

    init(tempPoints);
    tempPoints.kill();
}

MultipointEval::~MultipointEval() {
    clearSubProductTree(subProductTree);
}

void MultipointEval::initSubProductTree(Vec<Vec<ZZ_pEX>> &tree) {
    long k = NumBits(numPoints);
    long length = numPoints;
    tree.SetLength(k);

    for (long i = 0; i < k; i++) {
        tree[i].SetLength(length);
        length /= 2;
    }
}

void MultipointEval::clearSubProductTree(Vec<Vec<ZZ_pEX>> &tree) {
    for (long i = 0; i < tree.length(); i++)
        tree[i].kill();
    tree.kill();
}

void MultipointEval::buildSubProductTree(const Vec<ZZ_pE>& points) {
    long k = NumBits(numPoints);
    long length = numPoints;

    ZZ_pE one;
    ZZ_pE temp_coeff;

    set(one);

    // fill the first level as X - points[i]
    for (long i = 0; i < numPoints; i++) {
        NTL::negate(temp_coeff, points[i]);
        SetCoeff(subProductTree[0][i], 0, temp_coeff);
        SetCoeff(subProductTree[0][i], 1, one);
    }

    for (long i = 1; i < k; i++) {
        for (long j = 0; j < length / 2; j++) {
            mul(subProductTree[i][j], 
                subProductTree[i - 1][2 * j], 
                subProductTree[i - 1][2 * j + 1]);
        }

        if (length % 2 != 0)
            mul(subProductTree[i][length / 2 - 1], 
                subProductTree[i][length / 2 - 1],
                subProductTree[i - 1][length - 1]);

        length /= 2;
    }
}

void MultipointEval::goDownSubProductTree(Vec<ZZ_pE>& results,
                                          const ZZ_pEX& f) {
    long k = NumBits(numPoints);
    long length = 1;

    // reverse the bits of numPoints
    Util util;
    long revLength = util.reverseBits(numPoints);

    Vec<Vec<ZZ_pEX>> tempTree = subProductTree;
    rem(tempTree[k - 1][0], f, tempTree[k - 1][0]);
    revLength >>= 1;

    for (long j = k - 2; j >= 0; j--) {
        for (long i = 0; i < length; i++) {
            rem(tempTree[j][2 * i], 
                tempTree[j + 1][i], 
                tempTree[j][2 * i]);
            rem(tempTree[j][2 * i + 1], 
                tempTree[j + 1][i], 
                tempTree[j][2 * i + 1]);
        }

        if (revLength & 0x1)
            rem(tempTree[j][2 * length], 
                tempTree[j + 1][length - 1],
                tempTree[j][2 * length]);

        length = 2 * length + (revLength & 0x1);
        revLength >>= 1;
    }

    for (long i = 0; i < numPoints; i++)
        results[i] = coeff(tempTree[0][i], 0);
    
    clearSubProductTree(tempTree);
}

void MultipointEval::eval(Vec<ZZ_pE>& results,
                          const ZZ_pEX& f) {

    goDownSubProductTree(results, f);
}

void MultipointEval::eval(Vec<ZZ_pX>& results,
                          const ZZ_pEX& f) {

    Vec<ZZ_pE> tempResult;
    tempResult.SetLength(numPoints);

    eval(tempResult, f);

    for (long i = 0; i < numPoints; i++) {
        conv(results[i], tempResult[i]);
    }

    tempResult.kill();
}
