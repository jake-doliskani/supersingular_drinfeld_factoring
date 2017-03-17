
#ifndef FROBENIUSCOMP_H
#define FROBENIUSCOMP_H

#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>

using namespace NTL;

class FrobComp {
public:

    FrobComp();
    virtual ~FrobComp();

    /**
     * Computes x, tau^{-l}(x), tau^{-2l}(x), ..., tau^{-lm}(x) where tau is
     * the Frobenius x -> x^p.
     * @param result
     * @param m
     * @param l
     * @param initFrob the polynomial x^p mod {@code F}
     * @param F
     */
    void computeInverseFrobes(Vec<ZZ_pX>& result,
                              long m,
                              long l,
                              const ZZ_pX &initFrob,
                              const ZZ_pXModulus &F);

    /**
     * Given xqInit = x^{p^r}, this method computes x^{p^{2r}}, ..., x^{p^{mr}}.
     * @param result
     * @param m
     * @param xqInit
     * @param F
     */
    void computeFrobePowers(Vec<ZZ_pX>& result,
                            long m,
                            const ZZ_pX &xqInit,
                            const ZZ_pXModulus &F);

    /**
     * Computes f^{p^n} modulo {@code F}.
     * @param result
     * @param f
     * @param n
     * @param initFrob the polynomial x^p mod {@code F}
     * @param F
     */
    void computeFrobPower(ZZ_pX& result,
                          const ZZ_pX &f,
                          long n,
                          const ZZ_pX &initFrob,
                          const ZZ_pXModulus &F);


private:

};

#endif /* FROBENIUSCOMP_H */

