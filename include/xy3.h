#ifndef XY3_H
#define XY3_H


#include "clock3.h"


/* class : XY3
 * Class for handeling monte carlo simulations for the XY Model. Esentially the Clock2 class
 * with spins set to 50 (determined to converge to the XY Model at 50 spins).
 */
class XY3 : public Clock3
{
    public:
        XY3() = default;
        XY3(const int L);
        XY3(const XY3 &rhs);
        void set_spin();
        double sweep_energy(double beta, std::mt19937 &engine);
};


#endif
