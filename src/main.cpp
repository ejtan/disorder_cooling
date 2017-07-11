#include <iostream>
#include <vector>
#include <array>

#include "../include/neighbor.h"


/*-------------------------------------------------------------------------------------------------
 * GLOBAL CONSTANTS
 *-----------------------------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------------------------------
 * FUNCTIONS
 *-----------------------------------------------------------------------------------------------*/

/* test_neighbor_2D()
 * Tests sites for the 2D Neighbor table
 */
void test_neighbor_2D()
{
    std::cout << "  Testing 2D Neighbors... ";

    // Define array to hand code in the values expected for the neighbor indices.
    typedef std::array<int, 4> Neigh2;
    std::vector<Neigh2> expected_2D;

    // Push back expected neighbors, in order from position 0 to 8.
    // Values were hand computed.
    expected_2D.push_back( Neigh2{6, 1, 3, 2} ); // Row 0
    expected_2D.push_back( Neigh2{7, 2, 4, 0} );
    expected_2D.push_back( Neigh2{8, 0, 5, 1} );

    expected_2D.push_back( Neigh2{0, 4, 6, 5} ); // Row 1
    expected_2D.push_back( Neigh2{1, 5, 7, 3} );
    expected_2D.push_back( Neigh2{2, 3, 8, 4} );

    expected_2D.push_back( Neigh2{3, 7, 0, 8} ); // Row 2
    expected_2D.push_back( Neigh2{4, 8, 1, 6} );
    expected_2D.push_back( Neigh2{5, 6, 2, 7} );

    const int L = 3;
    bool isEqual = true;
    Neighbor<2> neigh;

    for (std::size_t i = 0; i < L * L; i++) {
        neigh.set_neighbors(i, L);

        if (neigh.neighbor != expected_2D[i]) {
            isEqual = false;
            break;
        } // Check if the generated and known index are true
    } // Loop over all lattice sites

    if (isEqual) std::cout << "Passed\n";
    else         std::cout << "Failed\n";
}


/* test_neighbor_3D()
 */
void test_neighbor_3D()
{
    std::cout << "  Testing 3D Neighbors... ";

    // Define array to hand code in the values expected for the neighbor indices.
    typedef std::array<int, 6> Neigh3;
    std::vector<Neigh3> expected_3D;

    // Push back expected neighbors in order from positons 0 to 26.
    // Values were hand computed.
    expected_3D.push_back( Neigh3{6, 1, 3, 2, 18, 9 } ); // Layer 0
    expected_3D.push_back( Neigh3{7, 2, 4, 0, 19, 10} );
    expected_3D.push_back( Neigh3{8, 0, 5, 1, 20, 11} );

    expected_3D.push_back( Neigh3{0, 4, 6, 5, 21, 12} );
    expected_3D.push_back( Neigh3{1, 5, 7, 3, 22, 13} );
    expected_3D.push_back( Neigh3{2, 3, 8, 4, 23, 14} );

    expected_3D.push_back( Neigh3{3, 7, 0, 8, 24, 15} );
    expected_3D.push_back( Neigh3{4, 8, 1, 6, 25, 16} );
    expected_3D.push_back( Neigh3{5, 6, 2, 7, 26, 17} );

    expected_3D.push_back( Neigh3{15, 10, 12, 11, 0, 18} ); // Layer 1
    expected_3D.push_back( Neigh3{16, 11, 13, 9 , 1, 19} );
    expected_3D.push_back( Neigh3{17, 9 , 14, 10, 2, 20} );

    expected_3D.push_back( Neigh3{9 , 13, 15, 14, 3, 21} );
    expected_3D.push_back( Neigh3{10, 14, 16, 12, 4, 22} );
    expected_3D.push_back( Neigh3{11, 12, 17, 13, 5, 23} );

    expected_3D.push_back( Neigh3{12, 16, 9 , 17, 6, 24} );
    expected_3D.push_back( Neigh3{13, 17, 10, 15, 7, 25} );
    expected_3D.push_back( Neigh3{14, 15, 11, 16, 8, 26} );

    expected_3D.push_back( Neigh3{24, 19, 21, 20, 9 , 0} ); // Layer 2
    expected_3D.push_back( Neigh3{25, 20, 22, 18, 10, 1} );
    expected_3D.push_back( Neigh3{26, 18, 23, 19, 11, 2} );

    expected_3D.push_back( Neigh3{18, 22, 24, 23, 12, 3} );
    expected_3D.push_back( Neigh3{19, 23, 25, 21, 13, 4} );
    expected_3D.push_back( Neigh3{20, 21, 26, 22, 14, 5} );

    expected_3D.push_back( Neigh3{21, 25, 18, 26, 15, 6} );
    expected_3D.push_back( Neigh3{22, 26, 19, 24, 16, 7} );
    expected_3D.push_back( Neigh3{23, 24, 20, 25, 17, 8} );

    const int L = 3;
    bool isEqual = true;
    Neighbor<3> neigh;

    for (std::size_t i = 0; i < L * L * L; i++) {
        neigh.set_neighbors(i, L);

        if (neigh.neighbor != expected_3D[i]) {
            isEqual = false;
            break;
        } // Check if the generated and known index are true
    } // Loop over all lattice sites

    if (isEqual) std::cout << "Passed\n";
    else         std::cout << "Failed\n";
}


/*-------------------------------------------------------------------------------------------------
 * MAIN
 *-----------------------------------------------------------------------------------------------*/

int main(int argc, char **argv)
{
    std::cout << "Testing Neighbor class implementation for correct neighbor indices\n";
    test_neighbor_2D();
    test_neighbor_3D();

    return 0;
}
