#ifndef XY2_H
#define XY2_H


#include "clock2.h"


/* class : XY2
 * Class for handeling monte carlo simulations for the XY Model. Esentially the Clock2 class
 * with spins set to 50 (determined to converge to the XY Model at 50 spins).
 */
class XY2 : public Clock2
{
    public:
        XY2() = default;
        XY2(const int L);
        XY2(const XY2 &rhs);
};


#endif
