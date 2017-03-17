
#include "BaseChange.h"

BaseChange::BaseChange() {
}

BaseChange::~BaseChange() {
}

void BaseChange::monomialToDual(Vec<ZZ_p>& dual,
                                const ZZ_pX& a,
                                const ZZ_pXModulus& modulus) {
    ZZ_pX temp;
    ZZ_pX modulusInvRev;

    long degree = deg(modulus);

    // compute 1 / rev(modulus, degree) mod x^degree
    reverse(modulusInvRev, modulus.f, degree);
    InvTrunc(modulusInvRev, modulusInvRev, degree);

    diff(temp, modulus);
    MulMod(temp, temp, a, modulus);
    reverse(temp, temp, degree - 1);
    MulTrunc(temp, temp, modulusInvRev, degree);

    for (long i = 0; i < degree; i++)
        GetCoeff(dual[i], temp, i);

    temp.kill();
    modulusInvRev.kill();
}

void BaseChange::dualToMonomial(ZZ_pX& result,
                                const Vec<ZZ_p>& dual,
                                const ZZ_pXModulus& modulus) {
    ZZ_pX temp1;
    ZZ_pX temp2;
    ZZ_pX modulusInv;

    long degree = deg(modulus);

    // compute 1 / modulus'
    diff(modulusInv, modulus.f);
    InvMod(modulusInv, modulusInv, modulus);

    for (long i = 0; i < degree; i++)
        SetCoeff(temp1, i, dual[i]);

    reverse(temp2, modulus.f, degree);
    MulTrunc(temp1, temp1, temp2, degree);
    reverse(temp2, temp1, degree - 1);

    MulMod(result, temp2, modulusInv, modulus);

    temp1.kill();
    temp2.kill();
    modulusInv.kill();
}

void BaseChange::changeBasis(ZZ_pX& result,
                             const ZZ_pX& f,
                             const ZZ_pX& g,
                             const ZZ_pXModulus& modulus) {

    ZZ_pX minPoly;
    long degree = deg(modulus);

    // compute the minimal polynomial of f
    MinPolyMod(minPoly, f, modulus);

    Vec<ZZ_p> dual;
    dual.SetLength(degree);
    clear(dual);

    // compute the dual basis image of x mod modulus
    monomialToDual(dual, g, modulus);

    // compute the power projection <dual, f>
    ProjectPowers(dual, dual, degree, f, modulus);

    ZZ_pXModulus tempModulus;
    build(tempModulus, minPoly);

    // compute result such that result(f) = g
    dualToMonomial(result, dual, tempModulus);

    minPoly.kill();
    dual.kill();
    tempModulus.f.kill();
}