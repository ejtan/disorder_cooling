#include "../include/model2.h"


/*-------------------------------------------------------------------------------------------------
 * PUBLIC METHODS
 *-----------------------------------------------------------------------------------------------*/

/* Default constructor
 */
Model2::Model2() : rand0(0.0, 1.0)
{
}


/* Constructor with parameters
 */
Model2::Model2(const int L) : rand0(0.0, 1.0)
{
    size = L * L;

    neigh.resize(size);
    J.resize(n_neigh * size);
}


/* Copy constructor
 */
Model2::Model2(const Model2 &rhs) : size(rhs.size), neigh(rhs.neigh), J(rhs.J)
{
}
