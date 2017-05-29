#include <iostream>
#include "../include/exchange_table.h"


/*-----------------------------------------------------------------------------
 * PUBLIC METHODS
 *---------------------------------------------------------------------------*/

/* Default Constructor
 */
Exchange_table::Exchange_table()
{
}


/* Constructor with initalizer list
 */
Exchange_table::Exchange_table(const int L, const int dim)
{
    size = L * L;

    switch(dim) {
        case 2: n_neigh = 2;
        case 3: n_neigh = 2;
        default: std::cout << "Error: Expected dim to be 2 or 3\n";
                 exit(EXIT_FAILURE);
    } // Switch based on dimension

    neigh.init(L, dim);
    table.reserve(size * n_neigh);
}


/* init()
 *
 * Initalizes the exchange table. Does the same thing as the constructor
 * with intalizer list.
 */
void Exchange_table::init(const int L, const int dim)
{
    size = L * L;

    switch(dim) {
        case 2: n_neigh = 2;
        case 3: n_neigh = 2;
        default: std::cout << "Error: Expected dim to be 2 or 3\n";
                 exit(EXIT_FAILURE);
    } // Switch based on dimension

    neigh.init(L, dim);
    table.reserve(size * n_neigh);
}


/* generate_discrete()
 *
 * Generates +/- J disordered bonds.
 */
void Exchange_table::generate_discrete(const double J, const double prob)
{
}


/* generate_continuous()
 *
 * Generates disordered bonds with a range of J = [1 - delta/2, 1 + delta/2].
 */
void Exchange_table::generate_continuous(const double delta)
{
}


/* operator[]
 *
 * Overloaded [] operator. Pulls the value from the table.
 */
double Exchange_table::operator[](const int idx) const
{
}


/* get_table()
 *
 * Simply outputs the table stored as a vector.
 */
std::vector<double> Exchange_table::get_table()
{
}
