#include <cmath>

#include "../include/clock.h"


/*------------------------------------------------------------------------------------------------
 * PRIVATE METHODS
 *----------------------------------------------------------------------------------------------*/

/* sweep_lattice()
 *
 * Performs Monte Carlo sweeps. Goes through the lattice once.
 * Has an equal chance of picking every site, so we randomlly pick a site L x L times.
 * Accepts / rejects flip based on the Metropolious Algorithm
 */
void Clock::sweep_lattice(const double beta, std::mt19937 &engine)
{
    int pos, new_angle;
    float delta_E;

    for (int i = 0; i < size; i++) {
        pos = static_cast<int>(rand0(engine) * size);

        // Compute new angle while the new angle is the same as the old one
        do
            new_angle = static_cast<int>(rand0(engine) * q);
        while (new_angle == spin[pos]);

        delta_E = calculate_deltaE(pos, new_angle);

        // Accept / reject
        if (rand0(engine) < exp(-beta * delta_E))
            spin[pos] = new_angle;
    } // Sweep over all possiable positions
}


/* calculate_deltaE()
 *
 * Computes the change in energy from changing the angle.
 */
double Clock::calculate_deltaE(const int pos, const int new_angle)
{
    int old_angle = spin[pos];
    int neigh_idx, neigh_angle, old_idx, new_idx;
    double E_tot = 0.0;

    for (int i = 0; i < n_neigh; i++) {
        // Get value of neighbor angle
        neigh_idx   = neigh[pos * n_neigh + i];
        neigh_angle = spin[neigh_idx];

        // Compute the index of old and new angle
        old_idx = (old_angle - neigh_angle + q) % q;
        new_idx = (new_angle - neigh_angle + q) % q;

        E_tot += J[i * n_neigh + i] * (cos_val[old_idx] - cos_val[new_idx]);
    } // Loop over neighboring spins

    // Return deltaE = E_old - E_new
    return E_tot;
}


/* calculate_totalE()
 *
 * Computes the total energy of the lattice.
 */
double Clock::calculate_totalE()
{
    double E = 0.0;
    int neigh_angle1, neigh_angle2, pos_angle;
    int E_idx1, E_idx2;

    for (int i = 0; i < size; i++) {
        neigh_angle1 = spin[neigh[i * n_neigh + 1]];
        neigh_angle2 = spin[neigh[i * n_neigh + 1]];
        pos_angle    = spin[i];

        E_idx1 = (pos_angle - neigh_angle1 + q) % q;
        E_idx2 = (pos_angle - neigh_angle2 + q) % q;

        E += -((J[i * n_neigh + 1] * cos_val[E_idx1]) +
               (J[i * n_neigh + 2] * cos_val[E_idx2]));
    } // Sweep over all sites

    return E;
}


/*------------------------------------------------------------------------------------------------
 * PUBLIC METHODS
 *----------------------------------------------------------------------------------------------*/

/* Default Constructor
 */
Clock::Clock() : rand0(0.0, 1.0)
{
}


/* Constructor with Initalizer list
 */
Clock::Clock(const int _L, const int _dim, const int _n_spin, const int _q) :
    Model(_L, _dim), rand0(0.0, 1.0), q(_q)
{
    spin.reserve(size);
    cos_val.reserve(_q);

    // Set spin table
    double dq = 2.0 * M_PI / static_cast<double>(_q);
    for (int i = 0; i < _q; i++)
        cos_val[i] = cos(i * dq);
}


/* init()
 *
 * Same as constructor with initalizer list
 */
void Clock::init(const int _L, const int _dim, const int _n_spin, const int _q)
{
    q = _q;

    Model::init(_L, _dim);
    spin.reserve(size);
    cos_val.reserve(_q);

    // Set spin table
    double dq = 2.0 * M_PI / static_cast<double>(_q);
    for (int i = 0; i < _q; i++)
        cos_val[i] = cos(i * dq);
}


/* set_spin()
 *
 * Sets the spin array with a value [0, q]. Index corresponds to the index of the
 * cos value table.
 */
void Clock::set_spin()
{
    std::random_device rd;
    std::mt19937 engine(rd());
    for (int i = 0; i < size; i++)
        spin[i] = static_cast<int>(rand0(engine) * q);
}


/* set_exchange()
 *
 * * Sets Exchange table with a mean of 1.
 */
void Clock::set_exchange(const double delta)
{
    J.generate_continuous(delta);
}


/* set_exchange()
 *
 * Sets exchange based on a discrete distribution.
 */
void Clock::set_exchange(const double J_val, const double p)
{
    J.generate_discrete(J_val, p);
}


/* sweep_energy()
 *
 * Performs Monte Carlo sweeps and measures the energy.
 */
double Clock::sweep_energy(const double beta, std::mt19937 &engine)
{
    double E_tot = 0.0;

    for (int i = 0; i < warmup; i++)
        sweep_lattice(beta, engine);

    for (int i = 0; i < measure; i++) {
        sweep_lattice(beta, engine);
        E_tot += calculate_totalE();
    } // Perform measurement sweeps

    return E_tot / static_cast<double>(measure * size);
}
