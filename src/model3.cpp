#include <random>

#include "../include/model3.h"


/*-------------------------------------------------------------------------------------------------
 * PUBLIC METHODS
 *-----------------------------------------------------------------------------------------------*/

/* Default constructor
 */
Model3::Model3() : isClean(true), rand0(0.0, 1.0)
{
}


/* Constructor with parameters
 */
Model3::Model3(const int L) : isClean(true), rand0(0.0, 1.0)
{
    size = L * L * L;

    neigh.resize(size);
    J.resize(n_neigh * size);

    for (size_t i = 0; i < size; i++)
        neigh[i].set_neighbors(i, L);
}


/* Copy constructor
 */
Model3::Model3(const Model3 &rhs) : size(rhs.size), isClean(rhs.isClean), neigh(rhs.neigh), J(rhs.J)
{
}


/* set_exchange()
 * Sets the exchange table used by 2D Models. delta is the range of the uniform distribution
 * with a mean centered at 1. The range of random values is J = [1 - delta/2, 1 + delta/2].
 */
void Model3::set_exchange(double delta)
{
    if (isClean)
        isClean = false;

    std::random_device rd;
    std::mt19937 engine(rd());
    double J_val, r_val;

    for (size_t i = 0; i < size; i++) {
        // 0 - 2 bond
        r_val = rand0(engine);
        if (rand0(engine) > 0.5) J_val = 1.0 - (delta * r_val / 2.0);
        else                     J_val = 1.0 + (delta * r_val / 2.0);
        J[i].J_arr[0]                    = J_val;
        J[neigh[i].neighbor[0]].J_arr[2] = J_val;

        // 1 - 3 bond
        r_val = rand0(engine);
        if (rand0(engine) > 0.5) J_val = 1.0 - (delta * r_val / 2.0);
        else                     J_val = 1.0 + (delta * r_val / 2.0);
        J[i].J_arr[1]                    = J_val;
        J[neigh[i].neighbor[1]].J_arr[3] = J_val;

        // 4 - 5 bond
        r_val = rand0(engine);
        if (rand0(engine) > 0.5) J_val = 1.0 - (delta * r_val / 2.0);
        else                     J_val = 1.0 + (delta * r_val / 2.0);
        J[i].J_arr[4]                    = J_val;
        J[neigh[i].neighbor[4]].J_arr[5] = J_val;
    } // Loop to set exchange table
}


/* get_exchange()
 * Returns the vector which contains the exchange table.
 */
std::vector<Exchange<3>> Model3::get_exchange() const
{
    return J;
}


/* set_run_param()
 * Overrides the default parameters for running simulations.
 */
void Model3::set_run_param(size_t Warmup, size_t Measure)
{
    warmup  = Warmup;
    measure = Measure;
}
