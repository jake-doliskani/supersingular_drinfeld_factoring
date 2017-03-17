
#ifndef BALANCEDMUL_H
#define BALANCEDMUL_H

#include <NTL/ZZ_pEX.h>
#include <NTL/matrix.h>

using namespace NTL;

class BalancedMul {
public:

    BalancedMul();
    virtual ~BalancedMul();

    /**
     * Computes the product of {@code vec[i]}, and put the result in 
     * {@code vec[0]}.
     * @param vec
     */
    void compute(Vec<ZZ_pEX> &vec);
    
    /**
     * Computes the product of {@code vec[i]} modulo {@code F}, and put the 
     * result in {@code vec[0]}.
     * @param vec
     */
    void compute(Vec<Mat<ZZ_pEX>>&vec, const ZZ_pEXModulus &F);

    /**
     * Computes the matrix product {@code result = A * B} modulo {@code F}.
     * @param result
     * @param A
     * @param B
     * @param F
     */
    void mul_ZZ_pXMat(Mat<ZZ_pX>& result,
                      const Mat<ZZ_pX>& A,
                      const Mat<ZZ_pX>& B,
                      const ZZ_pXModulus &F);

    /**
     * Computes the matrix product {@code result = A * B} modulo {@code F}.
     * @param result
     * @param A
     * @param B
     * @param F
     */
    void mul_ZZ_pEXMat(Mat<ZZ_pEX>& result,
                       const Mat<ZZ_pEX>& A,
                       const Mat<ZZ_pEX>& B,
                       const ZZ_pEXModulus &F);

    /**
     * Computes the matrix product B * A modulo {@code F}, where B has top row
     * [0 1] and bottom row {@code vec}.
     * @param result
     * @param vec
     * @param A
     * @param F
     */
    void mul_ZZ_pXMatSpec(Mat<ZZ_pX> &result,
                          const Vec<ZZ_pX> vec,
                          const Mat<ZZ_pX> A,
                          const ZZ_pXModulus &F);

    /**
     * Computes the matrix product B * A modulo {@code F}, where B has top row
     * [0 1] and bottom row {@code vec}.
     * @param result
     * @param vec
     * @param A
     * @param F
     */
    void mul_ZZ_pEXMatSpec(Mat<ZZ_pEX> &result,
                           const Vec<ZZ_pEX> vec,
                           const Mat<ZZ_pEX> A,
                           const ZZ_pEXModulus &F);

private:


};

#endif /* BALANCEDMUL_H */

