#ifndef EXCHANGE_H
#define EXCHANGE_H


#include <array>

/* struct : Exchange
 *
 * @template: Dim = Dimension of the lattice
 *
 * Provides bonds between neighboring spin sites. Follows convention from neighbor.h
 * Uses template to determine if 2D or 3D.
 *
 */
template <int Dim>
struct Exchange
{
    std::array<double, Dim * 2> J_arr;

    /* Constructors
     */
    Exchange() = default;
    Exchange(const Exchange<Dim> &rhs) : J_arr(rhs.J_arr) {}
};

#endif
