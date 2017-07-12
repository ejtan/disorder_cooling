#include "../include/model2.h"
#include <iostream>


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
