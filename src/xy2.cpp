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
