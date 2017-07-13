#include <random>

#include "../include/model2.h"


/*-------------------------------------------------------------------------------------------------
 * PUBLIC METHODS
 *-----------------------------------------------------------------------------------------------*/

/* Default constructor
 */
Model2::Model2() : rand0(0.0, 1.0), isClean(true)
{
}


/* Constructor with parameters
 */
Model2::Model2(const int L) : rand0(0.0, 1.0), isClean(true)
{
    size = L * L;

    neigh.resize(size);
    J.resize(n_neigh * size);

    for (size_t i = 0; i < size; i++)
        neigh[i].set_neighbors(i, L);
}


/* Copy constructor
 */
Model2::Model2(const Model2 &rhs) : size(rhs.size), isClean(rhs.isClean), neigh(rhs.neigh), J(rhs.J)
{
}


/* init()
 * Initalizes the model. Same as Constructor with parameters.
 */
void Model2::init(int L)
{
    size = L * L;

    neigh.resize(size);
    J.resize(n_neigh * size);

    for (size_t i = 0; i < size; i++)
        neigh[i].set_neighbors(i, L);
}


/* set_exchange()
 * Sets the exchange table used by 2D Models. delta is the range of the uniform distribution
 * with a mean centered at 1. The range of random values is J = [1 - delta/2, 1 + delta/2].
 */
void Model2::set_exchange(double delta)
{
    if (isClean)
        isClean = false;

    std::random_device rd;
    std::mt19937 engine(rd());
    double J_val, r_val;

    for (int i = 0; i < size; i++) {
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
    } // Loop to set exchange table
}


/* get_exchange()
 * Returns the vector which contains the exchange table.
 */
std::vector<Exchange<2>> Model2::get_exchange() const
{
    return J;
}
