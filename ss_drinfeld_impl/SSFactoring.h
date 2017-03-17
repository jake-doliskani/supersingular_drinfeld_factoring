#ifndef SSFACTORING_H
#define SSFACTORING_H


#include <NTL/ZZ_pXFactoring.h>


using namespace NTL;

class SSFactoring {
    
public:
    
    SSFactoring();
    virtual ~SSFactoring();
    
    /**
     * Factors the polynomial {@code f} using the proposed algorithm.
     * @param factors
     * @param f
     */
    void factor(Vec<pair_ZZ_pX_long> &factors, const ZZ_pX &f);
    
    /**
     * Recursively factors {@code f} using the {@link #split} method.
     * @param factors
     * @param f
     */
    void factorRec(Vec<ZZ_pX> &factors, const ZZ_pX &f);
    
    /**
     * Probabilistically computes a factor of {@code f} if it is reducible.
     * @param g a factor of {@code f}
     * @param f
     */
    void split(ZZ_pX &g, const ZZ_pX &f);
    
    /**
     * Sorts factors by multiplicity and degree.
     * @param factors
     */
    void sortFactors(Vec<pair_ZZ_pX_long> &factors);

private:

};

#endif /* SSFACTORING_H */

