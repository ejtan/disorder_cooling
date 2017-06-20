#include "../include/xy.h"


/*-------------------------------------------------------------------------------------------------
 * PRIVATE METHODS
 *-----------------------------------------------------------------------------------------------*/

/* Constructor with initalizer list
 */
XY::XY(const int L, const int dim) : Clock(L, dim, 50)
{
}


/* init()
 *
 * Same as constructor with initalizer list
 */
void XY::init(const int L, const int dim)
{
    Clock::init(L, dim, 50);
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
    Model::set_exchange(delta);
}


/* set_exchange()
 *
 * Sets discrete exchange.
 */
void XY::set_exchange(const double J_val, const double p)
{
    Model::set_exchange(J_val, p);
}


/* sweep_energy()
 *
 * Perform Monte Carlo sweeps and computes the energy.
 */
double XY::sweep_energy(const double beta, std::mt19937 &engine)
{
    Clock::sweep_energy(beta, engine);
}
