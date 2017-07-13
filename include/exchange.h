#ifndef EXCHANGE_H
#define EXCHANGE_H


#include <array>

/* struct : Exchange
 * Provides bonds between neighboring spin sites. Follows convention from neighbor.h
 * Uses template to determine if 2D or 3D.
 *
 */
template <std::size_t Dim>
struct Exchange
{
    std::array<double, Dim * 2> J_arr;

    /* Due to C++ compiler issues, implamentation is in header.
     */
    Exchange() = default;
    Exchange(const Exchange<Dim> &rhs) : J_arr(rhs.J_arr) {}
};

#endif
