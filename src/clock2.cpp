#include <cmath>
#include <random>

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


/* set_spin()
 * Sets the angle index representing the spin
 */
void Clock2::set_spin()
{
    std::random_device rd;
    std::mt19937 engine(rd());
    for (size_t i = 0; i < size; i++)
        spin[i] = static_cast<int>(rand0(engine) * q);
}
