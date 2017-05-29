#include <iostream>
#include "../include/neighbor_table.h"


/*-----------------------------------------------------------------------------
 * PRIVATE METHODS
 *---------------------------------------------------------------------------*/

/* set_2D_table()
 * Sets the neighbor table as a 2D square lattice.
 */
void set_2D_table(const int L)
{
}


/* set_2D_table()
 * Sets the neighbor table as a 3D square lattice.
 */
void set_3D_table(const int L)
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
    size = L * L;

    if (dim == 2) {
        n_neigh = 4;
        table.reserve(size * n_neigh);
        set_2D_table(L);
    } else if (dim == 3) {
        n_neigh = 6;
        table.reserve(size * n_neigh);
        set_3D_table(L);
    } else {
        std::cerr << "Error: Expected neighbor table dimension to be 2 or 3\n"
            << "Exiting Program." << std::endl;
        exit(EXIT_FAILURE);
    } // Perform check based on dimension input.
}


/* init()
 * Initalizes neighbor table like the constructor with initalizer list.
 */
void Neighbor_table::init(const int L, const int dim)
{
    size = L * L;

    if (dim == 2) {
        n_neigh = 4;
        table.reserve(size * n_neigh);
        set_2D_table(L);
    } else if (dim == 3) {
        n_neigh = 6;
        table.reserve(size * n_neigh);
        set_3D_table(L);
    } else {
        std::cerr << "Error: Expected neighbor table dimension to be 2 or 3\n"
            << "Exiting Program." << std::endl;
        exit(EXIT_FAILURE);
    } // Perform check based on dimension input.
}


/* operator[]
 * Returns the input index of the array.
 */
int Neighbor_table::operator[](const int idx)
{
}


/* get_table()
 * Outputs the neighbor table vector.
 */
std::vector<int> Neighbor_table::get_table()
{
}
