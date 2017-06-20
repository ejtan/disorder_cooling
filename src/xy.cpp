#include "../include/xy.h"


/*-------------------------------------------------------------------------------------------------
 * PRIVATE METHODS
 *-----------------------------------------------------------------------------------------------*/

/* Constructor
 */
XY::XY() : rand0(0.0, 1.0)
{
}


/* Constructor with initalizer list
 */
XY::XY(const int _L, const int _dim) : Clock(_L, _dim, 50, 50), rand0(0.0, 1.0)
{
}


/* init()
 *
 * Same as constructor with initalizer list
 */
void XY::init(const int _L, const int _dim)
{
    Clock::init(_L, _dim, 50, 50);
}


/* set_spin()
 *
 * Sets the spin of the XY Model
 */
void XY::set_spin()
{
    Clock::set_spin();
}


/* set_exchange()
 *
 * Sets continuous exchange centered at 1.
 */
void XY::set_exchange(const double delta)
{
    J.generate_continuous(delta);
}


/* set_exchange()
 *
 * Sets discrete exchange.
 */
void XY::set_exchange(const double J_val, const double p)
{
    J.generate_discrete(J_val, p);
}


/* sweep_energy()
 *
 * Perform Monte Carlo sweeps and computes the energy.
 */
double XY::sweep_energy(const double beta, std::mt19937 &engine)
{
    Clock::sweep_energy(beta, engine);
}
