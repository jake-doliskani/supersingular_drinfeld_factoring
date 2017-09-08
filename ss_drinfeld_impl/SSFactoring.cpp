

#include "SSFactoring.h"
#include "HasseLiftExt.h"
#include "Util.h"
#include <iostream>
#include <algorithm>

using namespace std;

SSFactoring::SSFactoring() {
}

SSFactoring::~SSFactoring() {
}

void SSFactoring::factor(Vec<pair_ZZ_pX_long>& factors, const ZZ_pX& f) {

    factors.SetLength(0);

    // separate the linear factors
    Vec<ZZ_pX> linearFactors;
    ZZ_pX temp;
    PowerXMod(temp, ZZ_p::modulus(), f);
    SetCoeff(temp, 1, coeff(temp, 1) - 1);
    GCD(temp, temp, f);
    RootEDF(linearFactors, temp);

    temp = f;
    for (long i = 0; i < linearFactors.length(); i++) {
        while (divide(temp, linearFactors[i])) {
            append(factors, cons(linearFactors[i], 1L));
            div(temp, temp, linearFactors[i]);
        }
    }

    cout << "linear factors: " << endl;
    cout << factors << endl;
    cout << "-----------------------" << endl;
    
    // now temp = f doesn't have any linear factors
    Vec<pair_ZZ_pX_long> sfd;
    SquareFreeDecomp(sfd, temp);

    Vec<ZZ_pX> tempFactors;
    for (long i = 0; i < sfd.length(); i++) {
        tempFactors.SetLength(0);
        factorRec(tempFactors, sfd[i].a);

        for (long j = 0; j < tempFactors.length(); j++)
            append(factors, cons(tempFactors[j], sfd[i].b));
    }

    linearFactors.kill();
    temp.kill();
}

void SSFactoring::sortFactors(Vec<pair_ZZ_pX_long>& factors) {
    Vec<pair_ZZ_pX_long> temp;
    temp.SetLength(0);
    long index = -1;

    for (long i = 0; i < factors.length(); i++) {
        for (long j = 0; j < temp.length(); j++) {
            if (factors[i] == temp[j]) {
                index = j;
                break;
            }
        }

        if (index >= 0)
            temp[index].b++;
        else
            temp.append(factors[i]);

        index = -1;
    }

    std::sort(temp.begin(), temp.end(),
              [](pair_ZZ_pX_long &a, pair_ZZ_pX_long & b) {
                  if (a.b > b.b)
                      return 0;
                  else if (a.b < b.b)
                      return 1;
                  if (deg(a.a) > deg(b.a))
                      return 0;
                  else
                      return 1;
              });

    factors = temp;
    temp.kill();
}

void SSFactoring::factorRec(Vec<ZZ_pX>& factors, const ZZ_pX& f) {
    if (ProbIrredTest(f)) {
        append(factors, f);
        return;
    }

    ZZ_pX temp;

    split(temp, f);
    while (deg(temp) == 0 || deg(temp) == deg(f))
        split(temp, f);

    cout << "f1: " << endl;
    cout << temp << endl;
    cout << "-----------------------" << endl;
    
    factorRec(factors, temp);
    div(temp, f, temp);

    cout << "f / f1: " << endl;
    cout << temp << endl;
    cout << "-----------------------" << endl;
    
    factorRec(factors, temp);

    temp.kill();
}

void SSFactoring::split(ZZ_pX& factor, const ZZ_pX& f) {
    ZZ_p a = random_ZZ_p();

    // dx = x - a
    ZZ_pX dx;
    SetX(dx);
    dx = dx - a;

    ZZ_pXModulus F;
    build(F, f);

    // temp = dx^{(p - 1) / 2} mod f
    ZZ_pX temp;
    PowerMod(temp, dx, (ZZ_p::modulus() - 1) / 2, F);

    // g = dx * (1 + dx^{(p - 1) / 2})
    ZZ_pX g = temp + 1;
    SqrMod(g, g, F);
    MulMod(g, g, dx, F);

    // delta = dx^{(p + 1) / 2} * (1 + dx^{(p - 1) / 2})^(p + 1)
    ZZ_pX delta = temp + 1;
    PowerMod(delta, delta, ZZ_p::modulus() + 1, F);
    MulMod(delta, delta, temp, F);
    MulMod(delta, delta, dx, F);

    // the Hasse lift
    Util util;
    ZZ_pX hasseLift;
    HasseLiftExt hasseLiftExt(g, delta, F);
    if (deg(f) < 700)
        hasseLiftExt.computeNaive(hasseLift, deg(f) / 2);
    else
        hasseLiftExt.compute(hasseLift, deg(f) / 2, 0);

    GCD(factor, hasseLift, f);

    if (deg(factor) > 0 && deg(factor) < deg(f)) {
        cout << "extension: " << endl;
        cout << dx << endl;
        cout << "-----------------------" << endl;

        cout << "elliptic module: " << endl;
        cout << "g:" << g << endl;
        cout << "delta:" << delta << endl;
        cout << "-----------------------" << endl;
    }

    dx.kill();
    temp.kill();
    g.kill();
    delta.kill();
    hasseLift.kill();
}

