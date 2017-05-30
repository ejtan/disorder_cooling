/* Written by : Eric Tan
 *
 * Test program to test neighbor table implamentation.
 * Compares output of a 2D 3x3 lattice with on computed by hand following
 * the convention set in include/neighbor_table.h
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include "../include/neighbor_table.h"


int main(void)
{
    Neighbor_table neigh(3, 2);
    std::vector<int> actual_2D = {6, 1, 3, 2, 7, 2, 4, 0, 8, 0, 5, 1, 0, 4, 6, 5,
        1, 5, 7, 4, 2, 3, 8, 4, 3, 7, 0, 8, 4, 8, 1, 6, 5, 6, 2, 7};
    std::vector<int> neigh_2D = neigh.get_table();

    if (neigh_2D.size() != actual_2D.size()) {
        std::cerr << "Error: Expected test and neighbor table to be same size"
            << std::endl;
        exit(EXIT_FAILURE);
    } // Check if both tables are the same size

    std::cout << "\nTesting if generated neighbor table is correct...\n\n";

    if (std::equal(actual_2D.begin(), actual_2D.end(), neigh_2D.begin()))
        std::cout << "2D Neighbor Table Test Passed.\n";
    else
        std::cout << "2D Neighbor Table Test Failed.\n";

    return 0;
}
