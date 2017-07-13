#include <cmath>
#include <boost/simd/exponential.hpp>
#include <boost/simd/function/multiplies.hpp>
#include "../include/ising2.h"


/*-------------------------------------------------------------------------------------------------
 * PRIVATE METHODS
 *-----------------------------------------------------------------------------------------------*/

/* sweep_lattice_clean()
 * Performans Monte Carlo sweeps. Sweeps the lattice once by choosing a random position and
 * proposing a spin flip using the Meteropolis Algorithm. This is done for the lattice size.
 */
void Ising2::sweep_lattice_clean(float beta, std::mt19937 &engine)
{
    namespace bs = boost::simd;

    for (size_t i = 0; i < size; i++) {
        int pos       = static_cast<int>(rand0(engine) * size);
        float delta_E = 2.0 * spin[pos] * (spin[neigh[pos].neighbor[0]] +
                                           spin[neigh[pos].neighbor[1]] +
                                           spin[neigh[pos].neighbor[2]] +
                                           spin[neigh[pos].neighbor[3]]);

        // Accept / reject flip
        // SIMD optimized functions used to compute exp(-beta * delta_E)
        if (rand0(engine) < bs::exp(bs::multiplies(-beta, delta_E)))
            spin[pos] = -spin[pos];
    } // Sweep over sites
}


/* sweep_lattice_clean()
 * Performans Monte Carlo sweeps. Sweeps the lattice once by choosing a random position and
 * proposing a spin flip using the Meteropolis Algorithm. This is done for the lattice size.
 */
void Ising2::sweep_lattice_disorder(float beta, std::mt19937 &engine)
{
    for (size_t i = 0; i < size; i++) {
        int pos       = static_cast<int>(rand0(engine) * size);
        float delta_E = 2.0 * spin[pos] * (J[pos].J_arr[0] * spin[neigh[pos].neighbor[0]] +
                                           J[pos].J_arr[1] * spin[neigh[pos].neighbor[1]] +
                                           J[pos].J_arr[2] * spin[neigh[pos].neighbor[2]] +
                                           J[pos].J_arr[3] * spin[neigh[pos].neighbor[3]]);

        // Accept / reject flip
        if (rand0(engine) < exp(-beta * delta_E))
            spin[pos] = -spin[pos];
    } // Sweep over sites
}


/*-------------------------------------------------------------------------------------------------
 * PUBLIC METHOD
 *-----------------------------------------------------------------------------------------------*/

/* Constructor with arguments
 */
Ising2::Ising2(const int L) : Model2(L)
{
    spin.resize(size);
}


/* Copy constructor
 */
Ising2::Ising2(const Ising2 &rhs) : Model2(rhs), spin(rhs.spin)
{
}


/* set_spin()
 * Sets the spin lattice to 1.
 */
void Ising2::set_spin()
{
    for (auto &&ele : spin)
        ele = 1;
}


/* sweep_energy()
 * Performs monte carlo sweeps and calcuates the energy
 */
double Ising2::sweep_energy(double beta, std::mt19937 &engine)
{
    double E_tot = 0.0;

    if (isClean) {
        for (size_t i = 0; i < warmup; i++)
            sweep_lattice_clean(beta, engine);

        for (size_t i = 0; i < measure; i++) {
            sweep_lattice_clean(beta, engine);

            #pragma omp simd reduction(+:E_tot)
            for (size_t j = 0; j < size; j++)
                E_tot += -spin[j] * (spin[neigh[j].neighbor[0]] + spin[neigh[j].neighbor[1]]);
        } // Perform measurement sweeps
    } else {
        for (size_t i = 0; i < warmup; i++)
            sweep_lattice_disorder(beta, engine);

        for (size_t i = 0; i < measure; i++) {
            sweep_lattice_disorder(beta, engine);

            for (size_t j = 0; j < size; j++)
                E_tot += -spin[j] * (J[j].J_arr[0] * spin[neigh[j].neighbor[0]] +
                                     J[j].J_arr[1] * spin[neigh[j].neighbor[1]]);
        } // Perform measurement sweeps
    } // Choose wheather to sweep with disorder or no disorder

    return E_tot / static_cast<double>(measure * size);
}
