#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <random>
#include <limits>

#include "../include/neighbor.h"
#include "../include/exchange.h"
#include "../include/ising2.h"
#include "../include/clock2.h"
#include "../include/xy2.h"
#include "../include/ising3.h"
#include "../include/clock3.h"
#include "../include/xy3.h"
#include "../include/disorder_cooling.h"


/*-------------------------------------------------------------------------------------------------
 * GLOBAL CONSTANTS
 *-----------------------------------------------------------------------------------------------*/
const int L     = 3;
const int N_pts = 200;
const int n_run = 10;
const double dT = 0.1;
const double delta = 5.0;


/*-------------------------------------------------------------------------------------------------
 * FORWARD DECLARATIONS
 *-----------------------------------------------------------------------------------------------*/
void test_neighbor_2D();
void test_neighbor_3D();
void test_exchange_2D();
void test_exchange_3D();
void test_ising(const std::array<double, N_pts> &T);
void test_clock(const std::array<double, N_pts> &T);
void test_xy(const std::array<double, N_pts> &T);


/*-------------------------------------------------------------------------------------------------
 * MAIN
 *-----------------------------------------------------------------------------------------------*/
int main(void)
{
    std::cout << "Testing Neighbor class implementation for correct neighbor indices\n";
    test_neighbor_2D();
    test_neighbor_3D();

    std::cout << "\nTesting for equality of exchange table\n";
    test_exchange_2D();
    std::cout << "\n";
    test_exchange_3D();

    // Initalize a temperature array
    std::array<double, N_pts> T;
    int curr = 0;

    std::generate(T.begin(), T.end(), [&curr]() {
            curr++;
            return dT * static_cast<double>(curr);
    });

    std::cout << "\nTesting 2D Ising Model with L = 4.\n";
    test_ising(T);

    std::cout << "\nTesting 2D Clock Model with L = 4.\n";
    test_clock(T);

    std::cout << "\nTesting 2D XY Mdoel with L = 4.\n";
    test_xy(T);

    return 0;
}


/*-------------------------------------------------------------------------------------------------
 * FUNCTIONS
 *-----------------------------------------------------------------------------------------------*/

/* test_exchange_2D()
 */
void test_exchange_2D()
{
    std::array<Neighbor<2>, L * L> neigh;
    Ising2 model(L);

    for (size_t i = 0; i < L * L; i++)
        neigh[i].set_neighbors(i, L);

    for (int i = 0; i < 5; i++) {
        std::cout << "  Testing 2D Bonds (test " << i + 1 << ")... ";

        std::random_device rd;
        std::mt19937 engine(rd());
        std::uniform_real_distribution<double> rand0(0.0, 20.0); // Choose a delta within this range

        model.set_exchange(rand0(engine));

        auto J = model.get_exchange();
        bool isEqual = true;

        for (int i = 0; i < L * L; i++) {
            double bond0 = J[i].J_arr[0];
            double bond2 = J[neigh[i].neighbor[0]].J_arr[2];
            double bond1 = J[i].J_arr[1];
            double bond3 = J[neigh[i].neighbor[1]].J_arr[3];

            if (fabs(bond0 - bond2) > std::numeric_limits<double>::epsilon() ||
                    fabs(bond1 - bond3) > std::numeric_limits<double>::epsilon()) {
                isEqual = false;
                break;
            } // Check if the bonds are the same
        } // Loop over sites

        if (isEqual) std::cout << "Passed\n";
        else         std::cout << "Failed\n";
    } // Loop over several tests
}


/* test_exchange_3D()
 */
void test_exchange_3D()
{
    std::array<Neighbor<3>, L * L> neigh;
    Ising3 model(L);

    for (size_t i = 0; i < L * L; i++)
        neigh[i].set_neighbors(i, L);

    for (int i = 0; i < 5; i++) {
        std::cout << "  Testing 3D Bonds (test " << i + 1 << ")... ";

        std::random_device rd;
        std::mt19937 engine(rd());
        std::uniform_real_distribution<double> rand0(0.0, 20.0); // Choose a delta within this range

        model.set_exchange(rand0(engine));

        auto J = model.get_exchange();
        bool isEqual = true;

        for (int i = 0; i < L * L; i++) {
            double bond0 = J[i].J_arr[0];
            double bond2 = J[neigh[i].neighbor[0]].J_arr[2];
            double bond1 = J[i].J_arr[1];
            double bond3 = J[neigh[i].neighbor[1]].J_arr[3];
            double bond4 = J[i].J_arr[4];
            double bond5 = J[neigh[i].neighbor[4]].J_arr[5];

            if (fabs(bond0 - bond2) > std::numeric_limits<double>::epsilon() ||
                    fabs(bond1 - bond3) > std::numeric_limits<double>::epsilon() ||
                    fabs(bond4 - bond5) > std::numeric_limits<double>::epsilon()) {
                isEqual = false;
                break;
            } // Check if the bonds are the same
        } // Loop over sites

        if (isEqual) std::cout << "Passed\n";
        else         std::cout << "Failed\n";
    } // Loop over several tests
}


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


/* test_ising()
 * Performs Monte carlo simulation for 2D clean system.
 */
void test_ising(const std::array<double, N_pts> &T)
{
    Ising2 ising2(4);
    Ising3 ising3(4);
    ising2.set_run_param(30000, 50000);
    ising3.set_run_param(30000, 50000);

    std::cout << "  Running 2D Ising Model (clean)... ";
    auto E_clean = compute_energy(T, ising2);
    compute_entropy(E_clean, T, 2, "2D_ising_clean.txt");
    write_energy(E_clean, T, "2D_ising_energy_clean.txt");
    std::cout << "Done\n";

    std::cout << "  Running 2D Ising Model (disorder)... ";
    auto E_disorder = compute_energy(T, ising2, delta, n_run);
    compute_entropy(E_disorder, T, 2, "2D_ising_disorder.txt");
    write_energy(E_disorder, T, "2D_ising_energy_disorder.txt");
    std::cout << "Done\n";

    std::cout << "  Running 3D Ising Model (clean)... ";
    E_clean = compute_energy(T, ising3);
    compute_entropy(E_clean, T, 2, "3D_ising_clean.txt");
    write_energy(E_clean, T, "3D_ising_energy_clean.txt");
    std::cout << "Done\n";

    std::cout << "  Running 3D Ising Model (disorder)... ";
    E_disorder = compute_energy(T, ising3, delta, n_run);
    compute_entropy(E_disorder, T, 2, "3D_ising_disorder.txt");
    write_energy(E_disorder, T, "3D_ising_energy_disorder.txt");
    std::cout << "Done\n";
}


/* test_clock()
 * Performs Monte carlo simulation for 2D clean system.
 */
void test_clock(const std::array<double, N_pts> &T)
{
    Clock2 clock2_2(4, 2), clock2_20(4, 20);
    Clock3 clock3_2(4, 2), clock3_20(4, 20);
    clock2_2.set_run_param(30000, 50000);
    clock2_20.set_run_param(30000, 50000);
    clock3_2.set_run_param(30000, 50000);
    clock3_20.set_run_param(30000, 50000);

    std::cout << "  Running 2D Clock model (clean) with 2 spins...  ";
    auto E2_clean = compute_energy(T, clock2_2);
    compute_entropy(E2_clean, T, 2, "2D_clock_clean_q=2.txt");
    write_energy(E2_clean, T, "2D_clock_energy_clean_q=2.txt");
    std::cout << "Done\n";

    std::cout << "  Running 2D Clock model (disorder) with 2 spins... ";
    auto E2_disorder = compute_energy(T, clock2_2, delta, n_run);
    compute_entropy(E2_disorder, T, 2, "2D_clock_disorder_q=2.txt");
    write_energy(E2_disorder, T, "2D_clock_energy_disorder_q=2.txt");
    std::cout << "Done\n";

    std::cout << "  Running 2D Clock model (clean) with 20 spins...  ";
    auto E20_clean = compute_energy(T, clock2_20);
    compute_entropy(E20_clean, T, 20, "2D_clock_clean_q=20.txt");
    write_energy(E20_clean, T, "2D_clock_energy_clean_q=20.txt");
    std::cout << "Done\n";

    std::cout << "  Running 2D Clock model (disorder) with 20 spins... ";
    auto E20_disorder = compute_energy(T, clock2_20, delta, n_run);
    compute_entropy(E20_disorder, T, 20, "2D_clock_disorder_q=20.txt");
    write_energy(E20_disorder, T, "2D_clock_energy_disorder_q=20.txt");
    std::cout << "Done\n";

    std::cout << "  Running 3D Clock model (clean) with 2 spins...  ";
    E2_clean = compute_energy(T, clock3_2);
    compute_entropy(E2_clean, T, 2, "3D_clock_clean_q=2.txt");
    write_energy(E2_clean, T, "3D_clock_energy_clean_q=2.txt");
    std::cout << "Done\n";

    std::cout << "  Running 3D Clock model (disorder) with 2 spins... ";
    E2_disorder = compute_energy(T, clock3_2, delta, n_run);
    compute_entropy(E2_disorder, T, 2, "3D_clock_disorder_q=2.txt");
    write_energy(E2_disorder, T, "3D_clock_energy_disorder_q=2.txt");
    std::cout << "Done\n";

    std::cout << "  Running 3D Clock model (clean) with 20 spins...  ";
    E20_clean = compute_energy(T, clock3_20);
    compute_entropy(E20_clean, T, 20, "3D_clock_clean_q=20.txt");
    write_energy(E20_clean, T, "3D_clock_energy_clean_q=20.txt");
    std::cout << "Done\n";

    std::cout << "  Running 3D Clock model (disorder) with 20 spins... ";
    E20_disorder = compute_energy(T, clock3_20, delta, n_run);
    compute_entropy(E20_disorder, T, 20, "3D_clock_disorder_q=20.txt");
    write_energy(E20_disorder, T, "3D_clock_energy_disorder_q=20.txt");
    std::cout << "Done\n";
}


/* test_xy()
 * Performs Monte carlo simulation for 2D clean system.
 */
void test_xy(const std::array<double, N_pts> &T)
{
    XY2 xy2(4);
    XY3 xy3(4);
    xy2.set_run_param(30000, 50000);
    xy3.set_run_param(30000, 50000);

    std::cout << "  Running 2D XY model (clean)... ";
    auto E_clean = compute_energy(T, xy2);
    compute_entropy(E_clean, T, 50, "2D_xy_clean.txt");
    write_energy(E_clean, T, "2D_xy_clean_energy.txt");
    std::cout << "Done\n";

    std::cout << "  Running 2D XY model (disorder)... ";
    auto E_disorder= compute_energy(T, xy2, delta, n_run);
    compute_entropy(E_disorder, T, 50, "2D_xy_disorder.txt");
    write_energy(E_disorder, T, "2D_xy_disorder_energy.txt");
    std::cout << "Done\n";

    std::cout << "  Running 3D XY model (clean)... ";
    E_clean = compute_energy(T, xy3);
    compute_entropy(E_clean, T, 50, "3D_xy_clean.txt");
    write_energy(E_clean, T, "3D_xy_clean_energy.txt");
    std::cout << "Done\n";

    std::cout << "  Running 3D XY model (disorder)... ";
    E_disorder= compute_energy(T, xy3, delta, n_run);
    compute_entropy(E_disorder, T, 50, "3D_xy_disorder.txt");
    write_energy(E_disorder, T, "3D_xy_disorder_energy.txt");
    std::cout << "Done\n";
}
