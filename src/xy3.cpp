#include "../include/xy3.h"


/*-------------------------------------------------------------------------------------------------
 * PUBLIC METHODS
 *-----------------------------------------------------------------------------------------------*/

/* Constructor with arguments
 */
XY3::XY3(const int L) : Clock3(L, 50)
{
}


/* Copy constructor
 */
XY3::XY3(const XY3 &rhs) : Clock3(rhs)
{
}


/* set_spin()
 * Calss clock model set_spin()
 */
void XY3::set_spin()
{
    Clock3::set_spin();
}


/* sweep_energy()
 * Calls clock mode sweep_energy()
 */
double XY3::sweep_energy(double beta, std::mt19937 &engine)
{
    return Clock3::sweep_energy(beta, engine);
}

/* sweep_binder()
 * Calls clock model sweep_binder()
 */
double XY3::sweep_binder(double beta, std::mt19937 &engine)
{
    return Clock3::sweep_binder(beta, engine);
}
