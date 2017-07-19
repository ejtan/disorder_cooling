#include <cmath>
#include <iostream>

#include "../include/clock2.h"


/*-------------------------------------------------------------------------------------------------
 * PUBLIC METHOD
 *-----------------------------------------------------------------------------------------------*/

/* sweep_lattice_clean()
 * Performans Monte Carlo sweeps. Sweeps the lattice once by choosing a random position and
 * proposing a spin flip using the Meteropolis Algorithm. This is done for the lattice size.
 */
void Clock2::sweep_lattice_clean(float beta, std::mt19937 &engine)
{
    for (size_t i = 0; i < size; i++) {
        size_t pos = static_cast<size_t>(rand0(engine) * size);

        // compute new angle
        int new_angle;
        do {
            new_angle = static_cast<int>(rand0(engine) * q);
        } while(new_angle == spin[pos]);

        // Compute the energy change
        double delta_E = 0.0;
        for (size_t i = 0; i < n_neigh; i++) {
            size_t neigh_angle = spin[neigh[pos].neighbor[i]];
            size_t old_angle   = spin[pos];

            size_t old_idx = (old_angle - neigh_angle + q) % q;
            size_t new_idx = (new_angle - neigh_angle + q) % q;

            delta_E += cos_val[old_idx] - cos_val[new_idx];
        } // Loop to compute total cos value

        if (rand0(engine) < exp(-beta * delta_E))
            spin[pos] = new_angle;
    } // Loop over sites
}


/* sweep_lattice_disorder()
 * Performans Monte Carlo sweeps. Sweeps the lattice once by choosing a random position and
 * proposing a spin flip using the Meteropolis Algorithm. This is done for the lattice size.
 */
void Clock2::sweep_lattice_disorder(float beta, std::mt19937 &engine)
{
    for (size_t i = 0; i < size; i++) {
        size_t pos = static_cast<size_t>(rand0(engine) * size);

        // Compute new angle
        int new_angle;
        do {
            new_angle = static_cast<int>(rand0(engine) * q);
        } while(new_angle == spin[pos]);

        // Compute energy change
        double delta_E = 0.0;
        for (size_t i = 0; i < n_neigh; i++) {
            size_t neigh_angle = spin[neigh[pos].neighbor[i]];
            size_t old_angle   = spin[pos];

            size_t old_idx = (old_angle - neigh_angle + q) % q;
            size_t new_idx = (new_angle - neigh_angle + q) % q;

            delta_E += J[pos].J_arr[i] * (cos_val[old_idx] - cos_val[new_idx]);
        }

        if (rand0(engine) < exp(-beta * delta_E))
            spin[pos] = new_angle;
    }
}


/*-------------------------------------------------------------------------------------------------
 * PUBLIC METHOD
 *-----------------------------------------------------------------------------------------------*/

/* Constructor with arguments
 */
Clock2::Clock2(const int L, const int _q) : Model2(L), q(_q)
{
    spin.resize(size);
    cos_val.resize(size);

    // Set spin table
    double dq = 2.0 * M_PI / static_cast<double>(_q);
    for (int i = 0; i < q; i++)
        cos_val[i] = cos(i * dq);
}


/* Copy constructor
 */
Clock2::Clock2(const Clock2 &rhs) : Model2(rhs),  q(rhs.q), spin(rhs.spin), cos_val(rhs.cos_val)
{
}


/* set_spin()
 * Sets the angle index representing the spin
 */
void Clock2::set_spin()
{
    std::random_device rd;
    std::mt19937 engine(rd());
    for (size_t i = 0; i < size; i++)
        spin[i] = static_cast<int>(rand0(engine) * q);
}


/* sweep_energy()
 */
double Clock2::sweep_energy(double beta, std::mt19937 &engine)
{
    double E_tot = 0.0;

    if (isClean) {
        for (size_t i = 0; i < warmup; i++)
            sweep_lattice_clean(beta, engine);

        for (size_t i = 0; i < measure; i++) {
            sweep_lattice_clean(beta, engine);

            for (size_t j = 0; j < size; j++) {
                size_t pos_angle = spin[j];
                size_t neigh1    = spin[neigh[j].neighbor[1]];
                size_t neigh2    = spin[neigh[j].neighbor[2]];

                size_t E_idx1 = (pos_angle - neigh1 + q) % q;
                size_t E_idx2 = (pos_angle - neigh2 + q) % q;

                E_tot += -(cos_val[E_idx1] + cos_val[E_idx2]);
            } // Compute energy of lattice
        } // Measurement sweeps
    } else {
        for (size_t i = 0; i < warmup; i++)
            sweep_lattice_disorder(beta, engine);

        for (size_t i = 0; i < measure; i++) {
            sweep_lattice_disorder(beta, engine);

            for (size_t j = 0; j < size; j++) {
                size_t pos_angle = spin[j];
                size_t neigh1    = spin[neigh[j].neighbor[1]];
                size_t neigh2    = spin[neigh[j].neighbor[2]];

                size_t E_idx1 = (pos_angle - neigh1 + q) % q;
                size_t E_idx2 = (pos_angle - neigh2 + q) % q;

                E_tot += -(J[j].J_arr[1] * cos_val[E_idx1] + J[j].J_arr[2] * cos_val[E_idx2]);
            } // Compute energy of lattice
        } // Measurement sweeps
    } // Choose if there is or isn't disorder

    return E_tot / static_cast<double>(measure * size);
}
