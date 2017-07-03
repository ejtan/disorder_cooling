#include <iostream>
#include "../include/neighbor_table.h"


/*-----------------------------------------------------------------------------
 * PRIVATE METHODS
 *---------------------------------------------------------------------------*/

/* set_2D_table()
 *
 * Sets the neighbor table as a 2D square lattice.
 */
void Neighbor_table::set_2D_table(const int L)
{
    for (int i = 0; i < N; i++) {
        // 0-th neighbor
        if (i - L < 0)
            table.push_back(i + N - L);
        else
            table.push_back(i - L);

        // 1-st neighbor
        if ((i + 1) % L == 0)
            table.push_back(i + 1 - L);
        else
            table.push_back(i + 1);

        // 2-nd neighbor
        if (i + L >= N)
            table.push_back(i + L - N);
        else
            table.push_back(i + L);

        // 3-rd neighbor
        if (i % L == 0)
            table.push_back(i + L - 1);
        else
            table.push_back(i - 1);
    } // Loop to set table
}


/* set_3D_table()
 *
 * Sets the neighbor table as a 3D square lattice.
 */
void Neighbor_table::set_3D_table(const int L)
{
}

/*-----------------------------------------------------------------------------
 * PUBLIC METHODS
 *---------------------------------------------------------------------------*/

/* Default constructor.
 */
Neighbor_table::Neighbor_table()
{
}


/* Constructor with initalizer list.
 */
Neighbor_table::Neighbor_table(const int L, const int dim)
{
    N = L * L;

    if (dim == 2) {
        n_neigh = 4;
        set_2D_table(L);
    } else if (dim == 3) {
        n_neigh = 6;
        set_3D_table(L);
    } else {
        std::cerr << "Error: Expected neighbor table dimension to be 2 or 3\n"
            << "Exiting Program." << std::endl;
        exit(EXIT_FAILURE);
    } // Perform check based on dimension input.
}


/* init()
 *
 * Initalizes neighbor table like the constructor with initalizer list.
 */
void Neighbor_table::init(const int L, const int dim)
{
    N = L * L;

    if (dim == 2) {
        n_neigh = 4;
        set_2D_table(L);
    } else if (dim == 3) {
        n_neigh = 6;
        set_3D_table(L);
    } else {
        std::cerr << "Error: Expected neighbor table dimension to be 2 or 3\n"
            << "Exiting Program." << std::endl;
        exit(EXIT_FAILURE);
    } // Perform check based on dimension input.
}


/* size()
 *
 * Returns the size of the neighbor table. This outputs the number of elements
 * stored in the table.
 * NOTE: For whatever reason, table.size() returns 0, so return N * n_neigh
 */
size_t Neighbor_table::size() const
{
    return N * n_neigh;
}


/* operator[]
 *
 * Returns the input index of the array.
 */
int Neighbor_table::operator[](const int idx) const
{
    return table[idx];
}


/* interator begin()
 */
std::vector<int>::iterator Neighbor_table::begin()
{
    return table.begin();
}


/* interator end()
 */
std::vector<int>::iterator Neighbor_table::end()
{
    return table.end();
}
