/* Written by : Eric Tan
 *
 * Test program to test neighbor table implamentation.
 * Compares output of a 2D 3x3 lattice with on computed by hand following
 * the convention set in include/neighbor_table.h
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include "../include/neighbor_table.h"
#include "../include/exchange_table.h"


/* test_neighbor()
 *
 * Compares output of a 2D 3x3 lattice with on computed by hand following
 * the convention set in include/neighbor_table.h
 */
void test_neighbor()
{
    Neighbor_table neigh(3, 2);

    // Hard coded hand computed neighbor table.
    std::vector<int> actual_2D = {6, 1, 3, 2, 7, 2, 4, 0, 8, 0, 5, 1, 0, 4, 6, 5,
        1, 5, 7, 3, 2, 3, 8, 4, 3, 7, 0, 8, 4, 8, 1, 6, 5, 6, 2, 7};

    if (neigh.size() != actual_2D.size()) {
        std::cerr << "Error: Expected test and neighbor table to be same size"
            << std::endl;
        exit(EXIT_FAILURE);
    } // Check if both tables are the same size

    std::cout << "Testing if generated neighbor table is correct... ";

    // Check if site table are equal
    if (std::equal(actual_2D.begin(), actual_2D.end(), neigh.begin()))
        std::cout << "2D Neighbor Table Test Passed.\n";
    else
        std::cout << "2D Neighbor Table Test Failed.\n";
}


/* check_exchange()
 *
 * Checks if connecting bonds are the same. Takes the difference and checks if
 * the absolute value is less than machine percision.
 */
bool check_exchange(const Exchange_table &J, const Neighbor_table &neigh)
{
    double bond1, bond2;

    for (int i = 0; i < 9; i++) {
        bond1 = J[i * 4];
        bond2 = J[neigh[i * 4] * 4 + 2];

        if (fabs(bond1 - bond2) > std::numeric_limits<double>::epsilon())
            return false;

        bond1 = J[i * 4 + 1];
        bond2 = J[neigh[i * 4 + 1] * 4 + 3];

        if (fabs(bond1 - bond2) > std::numeric_limits<double>::epsilon())
            return false;
    }

    return true;
}


/* test_exchange()
 *
 * Compares the bonds i - j to check if they are the same.
 */
void test_exchange()
{
    Exchange_table J(3, 2);
    Neighbor_table neigh(3, 2);

    std::cout << "Testing if continuous exchange table is correctly...\n";
    for (int i = 0; i < 5; i++) {
        std::cout << "  Performing Test " << i + 1 << "... ";

        J.generate_continuous(5.0);
        if (check_exchange(J, neigh))
            std::cout << "2D Exchange Table Test Passed.\n";
        else
            std::cout << "2D Exchange Table Test Failed.\n";
    }

    std::cout << "Testing if discrete exchange table is correctly...\n";
    for (int i = 0; i < 5; i++) {
        std::cout << "  Performing Test " << i + 1 << "... ";

        std::uniform_real_distribution<double> rand0(0.0, 1.0);
        std::random_device rd;
        std::mt19937 engine(rd());

        J.generate_discrete(1.0, rand0(engine));
        if (check_exchange(J, neigh))
            std::cout << "2D Exchange Table Test Passed.\n";
        else
            std::cout << "2D Exchange Table Test Failed.\n";
    }
}


int main(void)
{
    std::cout << "Performing Tests\n\n";
    test_neighbor();
    test_exchange();

    return 0;
}
