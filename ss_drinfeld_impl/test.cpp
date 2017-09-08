

#include <iostream>
#include <NTL/ZZ_pX.h>
#include <NTL/matrix.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZX.h>

#include "BaseChange.h"
#include "MultipointEval.h"
#include "Util.h"
#include "HasseLift.h"
#include "HasseLiftExt.h"
#include "MultiComposeMod.h"
#include "FrobComp.h"
#include "SSFactoring.h"



using namespace std;

//void testBaseChange(long degree) {
//    cout << "degree: " << degree << "\n";
//
//    ZZ_pX f, g, h;
//
//    random(f, degree - 1);
//    random(h, degree - 1);
//    
//    // build a square-free modulus
//    random(g, degree);
//    SetCoeff(g, degree, 1);
//    diff(h, g);
//    GCD(h, g, h);
//    div(g, g, h);
//    ZZ_pXModulus modulus;
//    build(modulus, g);
//    
//    CompMod(g, h, f, modulus);
//    BaseChange baseChane;
//    baseChane.changeBasis(f, f, g, modulus);
//
//    if (f == h)
//        cout << "ok\n";
//    else
//        cout << "oops\n";
//
//    f.kill();
//    g.kill();
//    h.kill();
//    modulus.f.kill();
//}
//

void testCompose() {
    long k = 1;
    long n = 2000;
    ZZ_pX f;
    random(f, n);
    SetCoeff(f, n, 1);
    ZZ_pXModulus F;
    build(F, f);

    Vec<ZZ_pX> result;
    Vec<ZZ_pX> result1;
    Vec<ZZ_pX> gVec;
    result.SetLength(k);
    result1.SetLength(k);
    gVec.SetLength(k);

    for (long i = 0; i < k; i++)
        random(gVec[i], n);

    ZZ_pX h;
    random(h, n);

    MultiComposeMod multiComposeMod;

    Util util;
    long start = util.getTimeMillis();
    multiComposeMod.compose(result, gVec, h, F);
    cout << util.getTimeMillis() - start << endl;

    start = util.getTimeMillis();
    for (long i = 0; i < k; i++) {
        CompMod(result1[i], gVec[i], h, F);
    }
    cout << util.getTimeMillis() - start << endl;

    for (long i = 0; i < k; i++) {
        if (result[i] != result1[i]) {
            cout << "Failed" << endl;
            break;
        }
    }
}
//
//
//void testMultipointEval() {
//    long numPoints = 100;
//    
//    long degree = 20;
//    ZZ_pX f;
//    random(f, degree);
//    SetCoeff(f, degree, 1);
//    
//    ZZ_pE::init(f);
//    
//    Vec<ZZ_pE> points;
//    Vec<ZZ_pE> results;
//    points.SetLength(numPoints);
//    results.SetLength(numPoints);
//    
//    for (long i = 0; i < numPoints; i++)
//        random(points[i]);
//    
//    ZZ_pEX g;
//    random(g, numPoints);
//    SetCoeff(g, numPoints + 1, 1);
//    
//    Util util;
//    long start = util.getTimeMillis();
//    MultipointEval multiPointEval;
//    multiPointEval.eval(results, g, points);
//    cout << (util.getTimeMillis() - start) << endl;
//    
//    start = util.getTimeMillis();
//    ZZ_pE temp;
//    for (long i = 0; i < numPoints; i++) {
//        eval(temp, g, points[i]);
//        if (temp != results[i]) {
//            cout << "not equal" << endl;
//            break;
//        }
//    }
//    cout << (util.getTimeMillis() - start) << endl;
//}
//

void testFrob() {
    ZZ_pX f;
    long d = 100;
    random(f, d);
    SetCoeff(f, d, 1);
    ZZ_pXModulus F;
    build(F, f);

    long m = 20;

    Vec<ZZ_pX> frobs;
    frobs.SetLength(m);
    ZZ_pX xq;
    SetX(xq);
    PowerMod(xq, xq, ZZ_p::modulus(), F);

    FrobComp frobComp;
    frobComp.computeFrobePowers(frobs, m, xq, F);

    Vec<ZZ_pX> frobs1;
    frobs1.SetLength(m);
    frobs1[0] = xq;
    for (long i = 1; i < m; i++) {
        PowerMod(frobs1[i], frobs1[i - 1], ZZ_p::modulus(), F);
    }

    if (frobs == frobs1)
        cout << "OK" << endl;
    else
        cout << "Failed" << endl;
}

void testHasseLift() {
    ZZ_pX f, g, delta, lift1, lift2;

    long degree = 10000;

    // build a square-free modulus
    Util util;
    util.randomMonic(f, degree);
    diff(g, f);
    GCD(g, g, f);
    div(f, f, g);
    ZZ_pXModulus F;
    build(F, f);
    degree = deg(F);

    random(g, degree);
    random(delta, degree);

    HasseLiftExt hasseLiftExt(g, delta, F);
    long start = util.getTimeMillis();
//    hasseLiftExt.computeNaive(lift1, degree / 2);
//    cout << util.getTimeMillis() - start << endl;
//
//    start = util.getTimeMillis();
    hasseLiftExt.compute(lift2, degree / 2, 1);
    cout << util.getTimeMillis() - start << endl;

    if (lift1 == lift2)
        cout << "OK" << endl;
    else
        cout << "Failed" << endl;
}

void testFactoring() {
    Util util;
    
    long n = 10;
    ZZ_pX f;
    util.randomMonic(f, n);
    cout << "p: " << ZZ_p::modulus() << endl;
    cout << "f: " << f << endl;
    
    Vec<pair_ZZ_pX_long> factors1;
    Vec<pair_ZZ_pX_long> factors2;
    
    SSFactoring sSfactoring;
//    long start = util.getTimeMillis();
    sSfactoring.factor(factors1, f);
    cout << factors1 << endl;
//    cout << util.getTimeMillis() - start << endl;
//    
//    start = util.getTimeMillis();
//    CanZass(factors2, f, 1);
//    cout << util.getTimeMillis() - start << endl;
//
//    ZZ_pX g;
//    mul(g, factors1);
//    if (f == g)
//        cout << "OK" << endl;
//    else
//        cout << "Failed" << endl;
}

int main(int argc, char** argv) {

    ZZ p;
    NextPrime(p, to_ZZ(7));
    ZZ_p::init(p);
    
    testFactoring();
    
    return 0;
}

