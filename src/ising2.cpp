#include "../include/ising2.h"


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
