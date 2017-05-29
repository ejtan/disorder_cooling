#include <iostream>
#include <random>
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
    if (prob > 1 || prob < 0) {
        std::cerr << "Error: Expected p to be between 0 and 1" << std::endl;
        exit(EXIT_FAILURE);
    } // Check if p is between 0 and 1

    std::uniform_real_distribution<float> rand0(0.0, 1.0);
    std::random_device rd;
    std::mt19937 engine(rd());
    double r, J_val;

    for (size_t i = 0; i < size; i++) {
        // 0 - 2 bond
        r = rand0(engine);
        if (r < prob) J_val = J;
        else          J_val = -J;
        table[i * n_neigh] = J_val;
        table[neigh[i * n_neigh] * n_neigh + 2] = J_val;

        // 1 - 3 bond
        r = rand0(engine);
        if (r < prob) J_val = J;
        else          J_val = -J;
        table[i * n_neigh + 1] = J_val;
        table[neigh[i * n_neigh + 1] * n_neigh + 3] = J_val;

        // 4 - 5 bond for 3D
        // Implament when 3D neighbor table is implamented
    } // Loop to set exchange table
}


/* generate_continuous()
 *
 * Generates disordered bonds with a range of J = [1 - delta/2, 1 + delta/2].
 */
void Exchange_table::generate_continuous(const double delta)
{
    std::uniform_real_distribution<float> rand0(0.0, 1.0);
    std::random_device rd;
    std::mt19937 engine(rd());
    double r, r_val , J_val;

    for (size_t i = 0; i < size; i++) {
        // 0 - 2 bond
        r = rand0(engine);
        r_val = rand0(engine);
        if (r > 0.5) J_val = 1.0 - (delta * r_val / 2.0);
        else         J_val = 1.0 - (delta * r_val / 2.0);
        table[i * n_neigh] = J_val;
        table[neigh[i * n_neigh] * n_neigh + 2] = J_val;

        // 1 - 3 bond
        r = rand0(engine);
        r_val = rand0(engine);
        if (r > 0.5) J_val = 1.0 - (delta * r_val / 2.0);
        else         J_val = 1.0 - (delta * r_val / 2.0);
        table[i * n_neigh + 1] = J_val;
        table[neigh[i * n_neigh + 1] * n_neigh + 3] = J_val;

        // 4 - 5 bond for 3D
        // Implament when 3D neighbor table is implamented
    } // Loop to set exchange table
}


/* operator[]
 *
 * Overloaded [] operator. Pulls the value from the table.
 */
double Exchange_table::operator[](const int idx) const
{
    return table[idx];
}


/* get_table()
 *
 * Simply outputs the table stored as a vector.
 */
std::vector<double> Exchange_table::get_table()
{
    return table;
}
