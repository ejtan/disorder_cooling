#include <cmath>

#include "../include/clock2.h"


/*-------------------------------------------------------------------------------------------------
 * PUBLIC METHOD
 *-----------------------------------------------------------------------------------------------*/

/* Constructor with arguments
 */
Clock2::Clock2(const int L, const int _q) : Model2(L), q(_q)
{
    spin.resize(size);
    cos_val.resize(size);

    // Set spin table
    double dq = 2.0 * M_PI / static_cast<double>(_q);
    for (int i = 0; i < q; i++)
        cos_val[i] = cos(i * dq);
}


/* Copy constructor
 */
Clock2::Clock2(const Clock2 &rhs) : Model2(rhs),  q(rhs.q), spin(rhs.spin), cos_val(rhs.cos_val)
{
}
