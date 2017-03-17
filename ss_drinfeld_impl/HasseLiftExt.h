
#ifndef HASSELIFT1_H
#define HASSELIFT1_H


#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vector.h>

using namespace NTL;

class HasseLiftExt {
public:

    /**
     * @param g first coefficient of the Drinfeld module.
     * @param delta second coefficient of the Drinfeld module
     * @param F the modulus
     */
    HasseLiftExt(const ZZ_pX &g,
                 const ZZ_pX &delta,
                 const ZZ_pXModulus &F);

    virtual ~HasseLiftExt();

    /**
     * Computes the n-th Hasse invariant lift using the proposed algorithm.
     * @param result
     * @param n
     * @param verbose
     */
    void compute(ZZ_pX &result, long n, int verbose = 0);

    /**
     * Computes the n-th Hasse invariant lift in a naive way.
     * @param result
     * @param n
     * @param verbose
     */
    void computeNaive(ZZ_pX& result, long n);

    
private:

    /**
     * Computes the matrix product tau^{l - 1}(A)tau^{l - 2}(A)...A.
     * @param result
     * @param A
     * @param l
     */
    void computeFrobProduct(Mat<ZZ_pEX> &result,
                            const Mat<ZZ_pEX> &A,
                            long l);

    /**
     * Computes the matrix product 
     * tau^{lm}(B(m))tau^{l(m - 1)}(B(m - 1))...B(0).
     * @param result
     * @param B
     * @param l
     */
    void computeFrobProduct(Mat<ZZ_pX> &result,
                            Mat<Vec<ZZ_pX>> &B,
                            long l);

    /**
     * Evaluates the matrix {@code A} at the polynomials {@code points}
     * @param result
     * @param A
     * @param points
     */
    void evaluate(Mat<Vec<ZZ_pX> > &result,
                  const Mat<ZZ_pEX> &A,
                  const Vec<ZZ_pX> &points);

    /**
     * Builds the bivariate polynomial needed for the multipoint evaluation.
     * @param inverseFrobs
     */
    void buildEvalModulus(const Vec<ZZ_pX> &inverseFrobs);

    /**
     * Builds the matrix A as proposed in the paper. The top row is [0 1] and
     * the bottom row is [delta * Y - tau(x) * delta, tau(g)].
     * @param result
     */
    void buildInitialMatrix(Mat<ZZ_pEX> &result);

    /**
     * Computes tau^{sl}(vec[s]), ..., tau^{l}(vec[1]), vec[0] where 
     * {@code vec.length() = s + 1}, and tau is the Frobenius x -> x^p.
     * @param vec
     * @param frobPowers
     */
    void actFrobPowers(Vec<ZZ_pX> &vec, const Vec<ZZ_pX> &frobPowers);

    /**
     * The baby-step-giant-step parameter
     */
    const double BETA = 0.8;
    
    /**
     * x^p
     */
    ZZ_pX frobenius;
    
    /**
     * The Drinfeld module coefficients
     */
    ZZ_pX g;
    ZZ_pX delta;
    
    
    ZZ_pXModulus F;
    ZZ_pEXModulus evalModulus;
};

#endif /* HASSELIFT_H */

