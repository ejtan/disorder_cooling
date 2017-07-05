#include <iostream>
#include <boost/simd/exponential.hpp>
#include <boost/simd/function/multiplies.hpp>
#include "../include/ising.h"


/*-------------------------------------------------------------------------------------------------
 * PRIVATE METHODS
 *-----------------------------------------------------------------------------------------------*/

/* sweep_lattice_disordered()
 *
 * Sweeps over the lattice of a disordered Ising system.
 */
void Ising::sweep_lattice_disordered(const float beta, std::mt19937 &engine)
{
    for (int i = 0; i < size; i++) {
        int pos       = static_cast<int>(rand0(engine) * size);
        float r_val   = rand0(engine);
        float delta_E = 0.0;

        for (int j = 0; j < n_neigh; j++)
            delta_E += J[pos * n_neigh + j] * spin[neigh[pos * n_neigh + j]];

        delta_E *= 2.0 * spin[pos];

        // Accept / reject flip
        // SIMD optimized with boost::simd (found to save a noticable amount of time).
        // Computes r_val < exp(-beta * delta_E)
        if (r_val < boost::simd::exp(boost::simd::multiplies(-beta, delta_E)))
            spin[pos] = -spin[pos];
    } // Sweep over all possiable sites
}


/* sweep_lattice_clean()
 *
 * Sweeps over a clean lattice system.
 */
void Ising::sweep_lattice_clean(const float beta, std::mt19937 &engine)
{
    for (int i = 0; i < size; i++) {
        int pos       = static_cast<int>(rand0(engine) * size);
        float r_val   = rand0(engine);
        float delta_E = 0.0;

        for (int j = 0; j < n_neigh; j++)
            delta_E += spin[neigh[pos * n_neigh + j]];

        delta_E *= 2.0 * spin[pos];

        // Accept / reject flip
        // SIMD optimized with boost::simd (found to save a noticable amount of time).
        // Computes r_val < exp(-beta * delta_E)
        if (r_val < boost::simd::exp(boost::simd::multiplies(-beta, delta_E)))
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

    if (J.is_clean()) {
        for (int i = 0; i < warmup; i++)
            sweep_lattice_clean(static_cast<float>(beta), engine);

        for (int i = 0; i < measure; i++) {
            sweep_lattice_clean(static_cast<float>(beta), engine);

            for (int j = 0; j < size; j++)
                E_tot += -spin[j] * (spin[neigh[j * n_neigh]] + spin[neigh[j * n_neigh + 1]]);
        } // Perform measurement sweeps
    } else {
        for (int i = 0; i < warmup; i++)
            sweep_lattice_disordered(static_cast<float>(beta), engine);

        for (int i = 0; i < measure; i++) {
            sweep_lattice_disordered(static_cast<float>(beta), engine);

            for (int j = 0; j < size; j++) {
                E_tot += -((J[j * n_neigh]     * spin[j] * spin[neigh[j * n_neigh]]) +
                           (J[j * n_neigh + 1] * spin[j] * spin[neigh[j * n_neigh + 1]]));
            } // Compute energy using 0 and 1 bonds.
        } // Perform measurement sweeps
    } // Function calls depending on if system is clean

    return E_tot / static_cast<double>(measure * size);
}
