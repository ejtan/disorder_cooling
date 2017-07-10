#ifndef NEIGHBOR_H
#define NEIGHBOR_H


#include <array>

/* struct : Neighbor
 * Provides neighboring indeces to a lattice site. Uses template to determine the dimension of the
 * lattice. Only supports 2D and 3D tables.
 * 2D:
 *      0
 *      |
 *  3 - x - 1
 *      |
 *      2
 *
 * 3D:
 *     4
 *     | 0
 *     |/
 *  3--x --1
 *    /|
 *   2 |
 *     5
 */
template <std::size_t Dim>
struct Neighbor
{
    static constexpr N = Dim * 2;
    std::array<int, N> neighbor;

    Neighbor();
    Neighbor(const Neighbor &rhs);
    set_neighbors(int pos);
};

#include "../src/neighbor.cpp"

#endif
