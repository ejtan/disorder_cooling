#include <iostream>
#include <vector>
#include <algorithm>

#include "models/ising.h"
#include "models/clock.h"
#include "models/xy.h"
#include "compute/energy.h"


/*-------------------------------------------------------------------------------------------------
 * GLOBAL CONSTANTS
 *-----------------------------------------------------------------------------------------------*/
const double delta = 5.0;
const int N_run = 10;


/*-------------------------------------------------------------------------------------------------
 * FORWARD DECLARATIONS
 *-----------------------------------------------------------------------------------------------*/
void test_ising(const std::vector<double> &T);
void test_clock2(const std::vector<double> &T);
void test_xy(const std::vector<double> &T);


/*-------------------------------------------------------------------------------------------------
 * MAIN
 *-----------------------------------------------------------------------------------------------*/
int main(void)
{
    const int N = 100;
    const double dT = 0.2;
    std::vector<double> T(N);

    std::generate(T.begin(), T.end(), [dT, curr=0]() mutable {
            curr++;
            return dT * static_cast<double>(curr);
    });

    std::cout << "\nTesting Ising.\n";
    test_ising(T);

    std::cout << "\nTesting Clock with 2 discrete spins.\n";
    test_clock2(T);

    std::cout << "\nTesting XY .\n";
    test_xy(T);

    std::cout << std::endl;
}


/*-------------------------------------------------------------------------------------------------
 * MAIN
 *-----------------------------------------------------------------------------------------------*/

/* test_ising()
 */
void test_ising(const std::vector<double> &T)
{
    // 2D Model
    Ising<2, 4> ising2;
    ising2.set_run_param(30000, 50000);

    std::cout << "\tRunning 2D Clean model..." << std::endl;
    auto E_clean2 = compute_energy(T, ising2);
    write_energy("ising2_clean_energy.txt", E_clean2, T);
    write_entropy("ising2_clean_entropy.txt", E_clean2, T, 2);

    std::cout << "\tRunning 2D Disorder model..." << std::endl;
    auto E_dis2 = compute_energy(T, ising2, delta, N_run);
    write_energy("ising2_disorder_energy.txt", E_dis2, T);
    write_entropy("ising2_disorder_entropy.txt", E_dis2, T, 2);

    // 3D Model
    Ising<3, 4> ising3;
    ising3.set_run_param(3000, 50000);

    std::cout << "\tRunning 3D Clean model..." << std::endl;
    auto E_clean3 = compute_energy(T, ising3);
    write_energy("ising3_clean_energy.txt", E_clean3, T);
    write_entropy("ising3_clean_entropy.txt", E_clean3, T, 2);

    std::cout << "\tRunning 3D Disorder model..." << std::endl;
    auto E_dis3 = compute_energy(T, ising3, delta, N_run);
    write_energy("ising3_disorder_energy.txt", E_dis3, T);
    write_entropy("ising3_disorder_entropy.txt", E_dis3, T, 2);
}


/* test_clock2()
 */
void test_clock2(const std::vector<double> &T)
{
    // 2D Model
    Clock<2, 4, 2> clock2;
    clock2.set_run_param(30000, 50000);

    std::cout << "\tRunning 2D clean model with 2 spins..." << std::endl;
    auto E_clean2 = compute_energy(T, clock2);
    write_energy("clock2_clean_energy.txt", E_clean2, T);
    write_entropy("clock2_clean_entropy.txt", E_clean2, T, 2);

    std::cout << "\tRunning 2D disorder model with 2 spins..." << std::endl;
    auto E_dis2 = compute_energy(T, clock2, delta, N_run);
    write_energy("clock2_disorder_energy.txt", E_dis2, T);
    write_entropy("clock2_disorder_entropy.txt", E_dis2, T, 2);

    // 3D Model
    Clock<3, 4, 2> clock3;
    clock2.set_run_param(30000, 50000);

    std::cout << "\tRunning 3D clean model with 2 spins..." << std::endl;
    auto E_clean3 = compute_energy(T, clock3);
    write_energy("clock3_clean_energy.txt", E_clean3, T);
    write_entropy("clock3_clean_entropy.txt", E_clean3, T, 2);

    std::cout << "\tRunning 3D disorder model with 2 spins..." << std::endl;
    auto E_dis3 = compute_energy(T, clock3, delta, N_run);
    write_energy("clock3_disorder_energy.txt", E_dis3, T);
    write_entropy("clock3_disorder_entropy.txt", E_dis3, T, 2);
}


/* test_xy()
 */
void test_xy(const std::vector<double> &T)
{
    // 2D Model
    XY<2, 4> xy2;
    xy2.set_run_param(30000, 50000);

    std::cout << "\tRunning 2D clean model..." << std::endl;
    auto E_clean2 = compute_energy(T, xy2);
    write_energy("xy2_clean_energy.txt", E_clean2, T);
    write_entropy("xy2_clean_entropy.txt", E_clean2, T, 50);

    std::cout << "\tRunning 2D disorder model..." << std::endl;
    auto E_dis2 = compute_energy(T, xy2, delta, N_run);
    write_energy("xy2_disorder_energy.txt", E_dis2, T);
    write_entropy("xy2_disorder_entropy.txt", E_dis2, T, 50);

    // 3D Model
    XY<3, 4> xy3;
    xy3.set_run_param(30000, 50000);

    std::cout << "\tRunning 3D clean model..." << std::endl;
    auto E_clean3 = compute_energy(T, xy3);
    write_energy("xy3_clean_energy.txt", E_clean3, T);
    write_entropy("xy3_clean_entropy.txt", E_clean3, T, 50);

    std::cout << "\tRunning 3D disorder model..." << std::endl;
    auto E_dis3 = compute_energy(T, xy3, delta, N_run);
    write_energy("xy3_disorder_energy.txt", E_dis3, T);
    write_entropy("xy3_disorder_entropy.txt", E_dis3, T, 50);

}
