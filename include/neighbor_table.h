#ifndef TABLE_H
#define TABLE_H

#include <vector>


/* Class : Neighbor_table
 *
 * Provides a neighbor table used to access neighbors of a given site in
 * a lattice. A square lattice is the only available in both 2 and 3 dimensions.
 * The table is set as a 1D vector with elements accessed by
 * x site * number of neighbors + i-th neighbor.
 *
 * Table follows the convention:
 *
 * 2D:
 *      0
 *      |
 *  3 - x - 1
 *      |
 *      2
 *
 * 3D: Do later.
 */
class Neighbor_table
{
    private:
        int size, n_neigh;
        std::vector<int> table;
        void set_2D_table(const int L);
        void set_3D_table(const int L);

    public:
        Neighbor_table();
        Neighbor_table(const int L, const int dim);
        void init(const int L, const int dim);
        int operator[](const int idx);
        std::vector<int> get_table();
};


#endif
