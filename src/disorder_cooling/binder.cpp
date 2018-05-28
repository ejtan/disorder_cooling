#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>

#include "compute/data_matrix.h"
#include "compute/binder.h"
#include "models/ising.h"
#include "models/clock.h"
#include "models/xy.h"


/*-------------------------------------------------------------------------------------------------
 * GLOBAL VARIABLES
 *-----------------------------------------------------------------------------------------------*/
const int N_pts = 16;
const int N_run = 5;
const int N_L = 4;
const int L1 = 4;
const int L2 = 6;
const int L3 = 8;
const double dT = 0.25;
const double delta = 1.0;


/*-------------------------------------------------------------------------------------------------
 * FORWARD DECLARATION
 *-----------------------------------------------------------------------------------------------*/
void test_ising_2D(const std::vector<double> &T);
void test_ising_3D(const std::vector<double> &T);
void test_clock_2D(const std::vector<double> &T);
void test_clock_3D(const std::vector<double> &T);
void test_xy_2D(const std::vector<double> &T);
void test_xy_3D(const std::vector<double> &T);


/*-------------------------------------------------------------------------------------------------
 * MAIN
 *-----------------------------------------------------------------------------------------------*/
int main(void)
{
    std::vector<double> T(N_pts);

    std::generate(T.begin(), T.end(), [curr=0]() mutable {
            curr++;
            return 1.0 + (dT * static_cast<double>(curr));
    });

    std::cout << "Testing Ising.\n";
    test_ising_2D(T);
    test_ising_3D(T);

    std::cout << "Testing Clock with 2 spins.\n";
    test_clock_2D(T);
    test_clock_3D(T);

    std::cout << "Testing XY.\n";
    test_xy_2D(T);
    test_xy_3D(T);
}


/*-------------------------------------------------------------------------------------------------
 * FUNCTIONS
 *-----------------------------------------------------------------------------------------------*/
void test_ising_2D(const std::vector<double> &T)
{
    std::vector<double> binder;
    Data_matrix<N_pts, N_L> data_clean, data_disorder;

    data_clean.insert_data(T);
    data_disorder.insert_data(T);

    Ising<2, L1> ising1;
    Ising<2, L2> ising2;
    Ising<2, L3> ising3;

    ising1.set_run_param(30000, 50000);
    ising2.set_run_param(30000, 50000);
    ising3.set_run_param(30000, 50000);

    // Clean binder
    std::cout << "\tPerforming 2D clean...\n" << std::flush;
    binder = compute_binder(T, ising1);
    data_clean.insert_data(binder);
    binder = compute_binder(T, ising2);
    data_clean.insert_data(binder);
    binder = compute_binder(T, ising3);
    data_clean.insert_data(binder);

    // Disorder binder
    std::cout << "\tPerforming 2D disorder...\n" << std::flush;
    binder = compute_binder(T, ising1, delta, N_run);
    data_disorder.insert_data(binder);
    binder = compute_binder(T, ising2, delta, N_run);
    data_disorder.insert_data(binder);
    binder = compute_binder(T, ising3, delta, N_run);
    data_disorder.insert_data(binder);

    std::ofstream of_clean("ising2_clean_binder.txt"), of_disorder("ising2_disorder_binder.txt");

    of_clean << data_clean;
    of_disorder << data_disorder;

    of_clean.close();
    of_disorder.close();
}


void test_ising_3D(const std::vector<double> &T)
{
    std::vector<double> binder;
    Data_matrix<N_pts, N_L> data_clean, data_disorder;

    data_clean.insert_data(T);
    data_disorder.insert_data(T);

    Ising<3, L1> ising1;
    Ising<3, L2> ising2;
    Ising<3, L3> ising3;

    ising1.set_run_param(30000, 50000);
    ising2.set_run_param(30000, 50000);
    ising3.set_run_param(30000, 50000);

    // Clean binder
    std::cout << "\tPerforming 3D clean...\n" << std::flush;
    binder = compute_binder(T, ising1);
    data_clean.insert_data(binder);
    binder = compute_binder(T, ising2);
    data_clean.insert_data(binder);
    binder = compute_binder(T, ising3);
    data_clean.insert_data(binder);

    // Disorder binder
    std::cout << "\tPerforming 3D disorder...\n" << std::flush;
    binder = compute_binder(T, ising1, delta, N_run);
    data_disorder.insert_data(binder);
    binder = compute_binder(T, ising2, delta, N_run);
    data_disorder.insert_data(binder);
    binder = compute_binder(T, ising3, delta, N_run);
    data_disorder.insert_data(binder);

    std::ofstream of_clean("ising3_clean_binder.txt"), of_disorder("ising3_disorder_binder.txt");

    of_clean << data_clean;
    of_disorder << data_disorder;

    of_clean.close();
    of_disorder.close();
}


void test_clock_2D(const std::vector<double> &T)
{
    std::vector<double> binder;
    Data_matrix<N_pts, N_L> data_clean, data_disorder;

    data_clean.insert_data(T);
    data_disorder.insert_data(T);

    Clock<2, L1, 2> clock1;
    Clock<2, L2, 2> clock2;
    Clock<2, L3, 2> clock3;

    clock1.set_run_param(30000, 50000);
    clock2.set_run_param(30000, 50000);
    clock3.set_run_param(30000, 50000);

    std::cout << "\tPerforming 2D clean...\n" << std::flush;
    binder = compute_binder(T, clock1);
    data_clean.insert_data(binder);
    binder = compute_binder(T, clock2);
    data_clean.insert_data(binder);
    binder = compute_binder(T, clock3);
    data_clean.insert_data(binder);

    std::cout << "\tPerforming 2D disorder...\n" << std::flush;
    binder = compute_binder(T, clock1, delta, N_run);
    data_disorder.insert_data(binder);
    binder = compute_binder(T, clock2, delta, N_run);
    data_disorder.insert_data(binder);
    binder = compute_binder(T, clock3, delta, N_run);
    data_disorder.insert_data(binder);

    std::ofstream of_clean("clock2_clean_binder.txt"), of_disorder("clock2_disorder_binder.txt");

    of_clean << data_clean;
    of_disorder << data_disorder;

    of_clean.close();
    of_disorder.close();
}


void test_clock_3D(const std::vector<double> &T)
{
    std::vector<double> binder;
    Data_matrix<N_pts, N_L> data_clean, data_disorder;

    data_clean.insert_data(T);
    data_disorder.insert_data(T);

    Clock<3, L1, 2> clock1;
    Clock<3, L2, 2> clock2;
    Clock<3, L3, 2> clock3;

    clock1.set_run_param(30000, 50000);
    clock2.set_run_param(30000, 50000);
    clock3.set_run_param(30000, 50000);

    std::cout << "\tPerforming 3D clean...\n" << std::flush;
    binder = compute_binder(T, clock1);
    data_clean.insert_data(binder);
    binder = compute_binder(T, clock2);
    data_clean.insert_data(binder);
    binder = compute_binder(T, clock3);
    data_clean.insert_data(binder);

    std::cout << "\tPerforming 3D disorder...\n" << std::flush;
    binder = compute_binder(T, clock1, delta, N_run);
    data_clean.insert_data(binder);
    binder = compute_binder(T, clock2, delta, N_run);
    data_clean.insert_data(binder);
    binder = compute_binder(T, clock3, delta, N_run);
    data_clean.insert_data(binder);

    std::ofstream of_clean("clock3_clean_binder.txt"), of_disorder("clock3_disorder_binder.txt");

    of_clean << data_clean;
    of_disorder << data_disorder;

    of_clean.close();
    of_disorder.close();
}


void test_xy_2D(const std::vector<double> &T)
{
    std::vector<double> binder;
    Data_matrix<N_pts, N_L> data_clean, data_disorder;

    data_clean.insert_data(T);
    data_disorder.insert_data(T);

    XY<2, L1> xy1;
    XY<2, L2> xy2;
    XY<2, L3> xy3;

    xy1.set_run_param(30000, 50000);
    xy2.set_run_param(30000, 50000);
    xy3.set_run_param(30000, 50000);

    std::cout << "\tPerforming 2D clean...\n" << std::flush;
    binder = compute_binder(T, xy1);
    data_clean.insert_data(binder);
    binder = compute_binder(T, xy2);
    data_clean.insert_data(binder);
    binder = compute_binder(T, xy3);
    data_clean.insert_data(binder);

    std::cout << "\tPerfoming 2D disorder...\n" << std::flush;
    binder = compute_binder(T, xy1, delta, N_run);
    data_disorder.insert_data(binder);
    binder = compute_binder(T, xy2, delta, N_run);
    data_disorder.insert_data(binder);
    binder = compute_binder(T, xy3, delta, N_run);
    data_disorder.insert_data(binder);

    std::ofstream of_clean("xy2_clean_binder.txt"), of_disorder("xy2_disorder_binder.txt");

    of_clean << data_clean;
    of_disorder << data_disorder;

    of_clean.close();
    of_disorder.close();
}


void test_xy_3D(const std::vector<double> &T)
{
    std::vector<double> binder;
    Data_matrix<N_pts, N_L> data_clean, data_disorder;

    data_clean.insert_data(T);
    data_disorder.insert_data(T);

    XY<3, L1> xy1;
    XY<3, L2> xy2;
    XY<3, L3> xy3;

    xy1.set_run_param(30000, 50000);
    xy2.set_run_param(30000, 50000);
    xy3.set_run_param(30000, 50000);

    std::cout << "\tPerforming 3D clean...\n" << std::flush;
    binder = compute_binder(T, xy1);
    data_clean.insert_data(binder);
    binder = compute_binder(T, xy2);
    data_clean.insert_data(binder);
    binder = compute_binder(T, xy3);
    data_clean.insert_data(binder);

    std::cout << "\tPerfoming 3D disorder...\n" << std::flush;
    binder = compute_binder(T, xy1, delta, N_run);
    data_disorder.insert_data(binder);
    binder = compute_binder(T, xy2, delta, N_run);
    data_disorder.insert_data(binder);
    binder = compute_binder(T, xy3, delta, N_run);
    data_disorder.insert_data(binder);

    std::ofstream of_clean("xy3_clean_binder.txt"), of_disorder("xy3_disorder_binder.txt");

    of_clean << data_clean;
    of_disorder << data_disorder;

    of_clean.close();
    of_disorder.close();
}
