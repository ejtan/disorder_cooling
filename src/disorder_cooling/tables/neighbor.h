#ifndef NEIGHBOR_H
#define NEIGHBOR_H


#include <array>

/* struct : Neighbor
 *
 * @template: Dim = number of dimension of the lattice
 * @template: L = number of sites in a single dimension
 *
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
template <int Dim, int L>
struct Neighbor
{
    std::array<int, Dim * 2> neighbor;

    /* Constructors
     */
    Neighbor() = default;
    Neighbor(const Neighbor<Dim, L> &rhs) : neighbor(rhs.neighbor) {}

    /* set_neighbors(pos)
     *
     * @param: pos = position on the lattice.
     *
     * Sets the neighbors to position pos.
     */
    void set_neighbors(int pos);
};

#include "neighbor.cpp"

#endif
