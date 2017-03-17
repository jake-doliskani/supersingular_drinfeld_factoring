
#ifndef MULTIPOINTEVAL_H
#define MULTIPOINTEVAL_H

#include <NTL/vector.h>
#include <NTL/ZZ_pEX.h>


using namespace NTL;

class MultipointEval {
public:

    
    MultipointEval(const Vec<ZZ_pE>& points);
    MultipointEval(const Vec<ZZ_pX>& points);
    virtual ~MultipointEval();

    /**
     * Evaluates {@code f} at the points given in the constructor of this class.
     * @param results
     * @param f
     */
    void eval(Vec<ZZ_pE> &results,
              const ZZ_pEX &f);

    /**
     * Evaluates {@code f} at the points given in the constructor of this class.
     * @param results
     * @param f
     */
    void eval(Vec<ZZ_pX> &results,
              const ZZ_pEX &f);

private:

    void init(const Vec<ZZ_pE>& points);
    
    /**
     * Allocates memory for the subproduct tree {@code tree}.
     * @param tree
     */
    void initSubProductTree(Vec<Vec<ZZ_pEX>> &tree);

    /**
     * Frees the memory of the subproduct tree {@code tree}.
     * @param tree
     */
    void clearSubProductTree(Vec<Vec<ZZ_pEX>> &tree);

    /**
     * Builds the subproduct tree {@code tree} using {@code points}.
     * @param points
     */
    void buildSubProductTree(const Vec<ZZ_pE> &points);
    
    /**
     * Computes the remainder of {@code f} divided by the leaves of the 
     * subproduct tree using the going down algorithm. See
     * "von zur Gathen, Gerhard, Modern Computer Algebra" for more details.
     * @param results
     * @param f
     */
    void goDownSubProductTree(Vec<ZZ_pE> &results, const ZZ_pEX &f);
    
    
    long numPoints;
    Vec<Vec<ZZ_pEX>> subProductTree;
};

#endif /* MULTIPOINTEVAL_H */

