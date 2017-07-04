#include <iostream>
#include "../include/model.h"


/*-------------------------------------------------------------------------------------------------
 * PUBLIC METHODS
 *-----------------------------------------------------------------------------------------------*/

/* Default Constructor
 */
Model::Model() : rand0(0.0, 1.0)
{
}


/* Constructor with arguments.
 */
Model::Model(const int L, const int dim) : rand0(0.0, 1.0)
{
    switch(dim) {
        case 2: n_neigh = 4; break;
        case 3: n_neigh = 6; break;
        default: std::cout << "Error: Expected dim to be 2 or 3" << std::endl;
                 exit(EXIT_FAILURE);
    } // Switch based on dimensionz

    size = L * L;
    neigh.init(L, dim);
    J.init(L, dim);
}


/* init()
 *
 * Behaves just as the constructor with arguments.
 */
void Model::init(const int L, const int dim)
{
    switch(dim) {
        case 2: n_neigh = 4; break;
        case 3: n_neigh = 6; break;
        default: std::cout << "Error: Expected dim to be 2 or 3" << std::endl;
                 exit(EXIT_FAILURE);
    } // Switch based on dimensionz

    size = L * L;
    neigh.init(L, dim);
    J.init(L, dim);
}


/* generate_exchange()
 *
 * Sets exchange table with a continuous distribution centered at 1.
 */
void Model::set_exchange(const double delta)
{
    J.generate_continuous(delta);
}


/* generate_exchange()
 *
 * Sets exchange table with a continuous distribution centered at 1.
 */
void Model::set_exchange(const double J_val, const double p)
{
    J.generate_discrete(J_val, p);
}
