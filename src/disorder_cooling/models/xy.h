#ifndef XY_H
#define XY_H

#include "clock.h"

// Set number of spins
constexpr int q = 50;


/* Class: XY
 *
 * @template Dim = dimension of system.
 * @template L = Number of sites per dimension
 *
 * Class for XY Model. Uses the clock model with spins set to 50.
 */
template <int Dim, int L>
class XY : public Clock<Dim, L, q>
{
    public:
        XY() = default;
        XY(const XY &rhs);
        void set_spin();
        double sweep_energy(double beta, std::mt19937 &engine);
        double sweep_binder(double beta, std::mt19937 &engine);
};

#include "xy.cpp"


#endif
