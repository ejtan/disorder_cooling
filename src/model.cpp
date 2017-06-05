#include <iostream>
#include "../include/model.h"


/* Default Constructor
 */
Model::Model()
{
}


/* Constructor with arguments.
 */
Model::Model(const int L, const int dim)
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
    J.generate_clean();
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
    J.generate_clean();
}
