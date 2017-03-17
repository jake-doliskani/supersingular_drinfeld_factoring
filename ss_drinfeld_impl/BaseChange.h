
#ifndef BASECHANGE_H
#define BASECHANGE_H

#include <NTL/ZZ_pX.h>

using namespace std;
using namespace NTL;

class BaseChange {
public:

    BaseChange();
    virtual ~BaseChange();

    /**
     * Computes the dual representation of {@code a} modulo {@code modulus}.
     * @param dual
     * @param a
     * @param modulus
     */
    void monomialToDual(Vec<ZZ_p> &dual,
                        const ZZ_pX &a,
                        const ZZ_pXModulus &modulus);

    /**
     * Computes the monomial representation of {@code dual} modulo 
     * {@code modulus}.
     * @param result
     * @param dual
     * @param modulus
     */
    void dualToMonomial(ZZ_pX &result,
                        const Vec<ZZ_p> &dual,
                        const ZZ_pXModulus &modulus);

    /**
     * Given {@code f} and {@code g} polynomials in $\mathbb{F}_p[X]/(modulus)$,
     * computes a polynomial $h \in \mathbb{F}_p[X]/(modulus)$ such that
     * $h(f) = g$ if such polynomial exists.  
     * @param result	a polynomial h such that $h(f) = g$ if such 
     * polynomial exists
     */
    void changeBasis(ZZ_pX &result,
                     const ZZ_pX &f,
                     const ZZ_pX &g,
                     const ZZ_pXModulus &modulus);
private:


};

#endif /* BASECHANGE_H */

