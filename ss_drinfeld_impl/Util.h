#ifndef UTIL_H
#define UTIL_H

#include <NTL/ZZ_pX.h>

using namespace NTL;


class Util {
    
public:
    
    Util();
    virtual ~Util();

    /**
     * @return the system time in milliseconds
     */
    long getTimeMillis();
    
    /**
     * @param a
     * @return a with the bits in reverse order
     */
    long reverseBits(long a);
    
    /**
     * @param x a random monic polynomial of degree {@code degree}
     * @param degree
     */
    void randomMonic(ZZ_pX& x, int degree);

private:

};

#endif /* UTIL_H */

