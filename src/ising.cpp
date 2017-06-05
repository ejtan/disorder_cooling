#include <iostream>
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
    spin.reserve(size);
}


/* init()
 *
 * Same as constructor with initalizer list
 */
void Ising::init(const int L, const int dim)
{
    Model::init(L, dim);
    spin.reserve(size);
}


/* set_spin()
 *
 * Sets the spin of lattice. Initalizes to +1.
 */
void Ising::set_spin()
{
    if (size) {
        for (int i = 0; i < size; i++)
            spin[i] = 1;
    } else {
        std::cerr << "Ising class is not initalized" << std::endl;
        exit(EXIT_FAILURE);
    } // Check if Class is initalized
}


/* sweep_energy()
 *
 * Performs Monte Carlo sweeps and measures the energy
 */
double Ising::sweep_energy()
{
}
