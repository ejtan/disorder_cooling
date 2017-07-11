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
    std::array<int, Dim * 2> neighbor;
    
    /* Due to the way C++ handles compiling template structs/classes, the implamentation is
     * here instead of two seperate files (Check the best way to handle this later).
     */
    Neighbor() = default;

    Neighbor(const Neighbor<Dim> &rhs) : neighbor(rhs.neighbor) {}

    void set_neighbors(int pos, int L)
    {
        /* In C++17 standard, replace if statements with if constexpr(...)
         */
        if (Dim == 2 || Dim == 3) {
            int N2 = L * L; // Number of sites on a 2D layer

            // 0-th neighbor
            if ((pos % N2) - L < 0)   neighbor[0] = pos + N2 - L;
            else                      neighbor[0] = pos - L;
            // 1-st neighbor
            if ((pos + 1) % L == 0)   neighbor[1] = pos + 1 - L;
            else                      neighbor[1] = pos + 1;
            // 2-nd neighbor
            if ((pos + L) % N2 < L) neighbor[2] = pos + L - N2;
            else                      neighbor[2] = pos + L;
            // 3-rd neighbor
            if (pos % L == 0)         neighbor[3] = pos + L - 1;
            else                      neighbor[3] = pos - 1;

            if (Dim == 3) {
                int N3 = L * L * L; // Total number of sites in a 3D lattice

                // 4-th neighbor
                if (pos - N2 < 0)   neighbor[4] = pos + N3 - N2;
                else                neighbor[4] = pos - N2;
                // 5-th neighbor
                if (pos + N2 >= N3) neighbor[5] = pos - N3 + N2;
                else                neighbor[5] = pos + N2;
            } // Set 3D layer (if needed)
        } // Set 2D layer
    }
};

#endif
