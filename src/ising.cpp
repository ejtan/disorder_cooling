#include <iostream>
#include <boost/simd/exponential.hpp>
#include "../include/ising.h"


/*-------------------------------------------------------------------------------------------------
 * PRIVATE METHODS
 *-----------------------------------------------------------------------------------------------*/

/* sweep_lattice()
 *
 * Sweeps over the lattice of a disordered Ising system
 */
void Ising::sweep_lattice(const double beta, std::mt19937 &engine)
{
    int pos;
    float r_val, delta_E;

    for (int i = 0; i < size; i++) {
        pos     = static_cast<int>(rand0(engine) * size);
        r_val   = rand0(engine);
        delta_E = 2.0 * spin[pos] *
            ((J[pos * n_neigh]     * spin[neigh[pos * n_neigh]]) +
             (J[pos * n_neigh + 1] * spin[neigh[pos * n_neigh + 1]]) +
             (J[pos * n_neigh + 2] * spin[neigh[pos * n_neigh + 2]]) +
             (J[pos * n_neigh + 3] * spin[neigh[pos * n_neigh + 3]]));

        // Accept / reject flip
        if (r_val < boost::simd::exp(-beta * delta_E))
            spin[pos] = -spin[pos];
    } // Sweep over all possiable sites
}


/*-------------------------------------------------------------------------------------------------
 * PUBLIC METHODS
 *-----------------------------------------------------------------------------------------------*/

/* Constructor with Initalizer list
 */
Ising::Ising(const int L, const int dim) : Model(L, dim)
{
    spin.resize(size);
}


/* copy constructor
 */
Ising::Ising(const Ising &rhs) : Model(rhs), spin(rhs.spin)
{
}

/* init()
 *
 * Same as constructor with initalizer list
 */
void Ising::init(const int L, const int dim)
{
    Model::init(L, dim);
    spin.resize(size);
}


/* set_spin()
 *
 * Sets the spin of lattice. Initalizes to +1.
 */
void Ising::set_spin()
{
    for (int i = 0; i < size; i++)
        spin[i] = 1;
}


/* set_exchange()
 *
 * Sets exchange table with mean of 1.
 */
void Ising::set_exchange(const double delta)
{
    Model::set_exchange(delta);
}


/* set_exchange()
 *
 * Sets exchange table with discrete probability
 */
void Ising::set_exchange(const double J_val, const double p)
{
    Model::set_exchange(J_val, p);
}


/* sweep_energy()
 *
 * Performs Monte Carlo sweeps and measures the energy
 */
double Ising::sweep_energy(const double beta, std::mt19937 &engine)
{
    double E_tot = 0.0;

    for (int i = 0; i < warmup; i++)
        sweep_lattice(beta, engine);

    for (int i = 0; i < measure; i++) {
        sweep_lattice(beta, engine);

        for (int j = 0; j < size; j++) {
            E_tot += -((J[j * n_neigh]     * spin[j] * spin[neigh[j * n_neigh]]) +
                       (J[j * n_neigh + 1] * spin[j] * spin[neigh[j * n_neigh + 1]]));
        } // Compute energy using 0 and 1 bonds.
    } // Perform measurement sweeps

    return E_tot / static_cast<double>(measure * size);
}
