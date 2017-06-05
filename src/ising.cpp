#include "../include/ising.h"


/* Default Constructor
 */
Ising::Ising() : rand0(0.0, 1.0)
{
}


/* Constructor with Initalizer list
 */
Ising::Ising(const int L, const int dim) : Model(L, dim), rand0(0.0, 1.0)
{
}


/* init()
 *
 * Same as constructor with initalizer list
 */
void Ising::init(const int L, const int dim)
{
}


/* set_spin()
 *
 * Sets the spin of lattice. Initalizes to +/- 1.
 */
void Ising::set_spin()
{
}


/* sweep_energy()
 *
 * Performs Monte Carlo sweeps and measures the energy
 */
double Ising::sweep_energy()
{
}
