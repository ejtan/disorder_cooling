#include "../include/xy2.h"


/*-------------------------------------------------------------------------------------------------
 * PUBLIC METHODS
 *-----------------------------------------------------------------------------------------------*/

/* Constructor with arguments
 */
XY2::XY2(const int L) : Clock2(L, 50)
{
}


/* Copy constructor
 */
XY2::XY2(const XY2 &rhs) : Clock2(rhs)
{
}


/* set_spin()
 * Calss clock model set_spin()
 */
void XY2::set_spin()
{
    Clock2::set_spin();
}


/* sweep_energy()
 * Calls clock mode sweep_energy()
 */
double XY2::sweep_energy(double beta, std::mt19937 &engine)
{
    return Clock2::sweep_energy(beta, engine);
}
