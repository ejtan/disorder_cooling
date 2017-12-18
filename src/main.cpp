#include <iostream>
#include <array>
#include <algorithm>

#include "../include/disorder_cooling.h"


/*-------------------------------------------------------------------------------------------------
 * GLOBAL VARIABLES
 *-----------------------------------------------------------------------------------------------*/
const int N_L = 3;
const int N_pts = 60;
const int N_run = 3;
const std::array<int, N_L> L = {4, 6, 8};
const double dT = 0.1;
const double delta = 0.5;


/*-------------------------------------------------------------------------------------------------
 * FORWARD DECLARATIONS
 *-----------------------------------------------------------------------------------------------*/
void test_ising(const std::array<double, N_pts> &T);
void test_clock(const std::array<double, N_pts> &T);
void test_xy(const std::array<double, N_pts> &T);


/*-------------------------------------------------------------------------------------------------
 * MAIN
 *-----------------------------------------------------------------------------------------------*/
int main(void)
{
    std::array<double, N_pts> T;
    int curr = 0;

    std::generate(T.begin(), T.end(), [&curr]() {
            curr++;
            return dT * static_cast<double>(curr);
    });

    test_ising(T);
    test_clock(T);
    test_xy(T);
}


/*-------------------------------------------------------------------------------------------------
 * FUNCTIONS
 *-----------------------------------------------------------------------------------------------*/

/* test_ising()
 * Tests Ising model.
 */
void test_ising(const std::array<double, N_pts> &T)
{
    Data_matrix data_clean2(N_pts, N_L+1);
    Data_matrix data_clean3(N_pts, N_L+1);
    Data_matrix data_disorder2(N_pts, N_L+1);
    Data_matrix data_disorder3(N_pts, N_L+1);

    data_clean2.insert_array(T.data());
    data_clean3.insert_array(T.data());
    data_disorder2.insert_array(T.data());
    data_disorder3.insert_array(T.data());

    std::cout << "Performing 2D Ising (clean)\n";
    for (int i = 0; i < N_L; i++) {
        std::cout << "\tL = " << L[i] << "... ";
        Ising2 ising(L[i]);
        ising.set_run_param(30000, 50000);
        auto binder = compute_binder(T, ising);
        data_clean2.insert_array(binder.data());
    }

    std::cout << "Performing 2D Ising (disorder)\n";
    for (int i = 0; i < N_L; i++) {
        std::cout << "\tL = " << L[i] << "... ";
        Ising2 ising(L[i]);
        ising.set_run_param(30000, 50000);
        auto binder = compute_binder(T, ising, delta, N_run);
        data_disorder2.insert_array(binder.data());
    }

    std::cout << "Performing 3D Ising (clean)\n";
    for (int i = 0; i < N_L; i++) {
        std::cout << "\tL = " << L[i] << "... ";
        Ising3 ising(L[i]);
        ising.set_run_param(30000, 50000);
        auto binder = compute_binder(T, ising);
        data_clean3.insert_array(binder.data());
    }

    std::cout << "Performing 3D Ising (clean)\n";
    for (int i = 0; i < N_L; i++) {
        std::cout << "\tL = " << L[i] << "... ";
        Ising3 ising(L[i]);
        ising.set_run_param(30000, 50000);
        auto binder = compute_binder(T, ising, delta, N_run);
        data_disorder3.insert_array(binder.data());
    }

    std::ofstream of_clean2("binder_clean_ising2.dat"), of_disorder2("binder_disorder_ising2.dat"),
        of_clean3("binder_clean_ising3.dat"), of_disorder3("binder_disorder_ising3.dat");

    of_clean2 << data_clean2;
    of_clean3 << data_clean3;
    of_disorder2 << data_disorder2;
    of_disorder3 << data_disorder3;

    of_clean2.close();
    of_clean3.close();
    of_disorder2.close();
    of_disorder3.close();
}


/* test_clock()
 * Tests Clock Model.
 */
void test_clock(const std::array<double, N_pts> &T)
{
    Data_matrix data_clean2(N_pts, N_L+1);
    Data_matrix data_clean3(N_pts, N_L+1);
    Data_matrix data_disorder2(N_pts, N_L+1);
    Data_matrix data_disorder3(N_pts, N_L+1);

    data_clean2.insert_array(T.data());
    data_clean3.insert_array(T.data());
    data_disorder2.insert_array(T.data());
    data_disorder3.insert_array(T.data());

    std::cout << "Performing 2D clock (2 spins) (clean)\n";
    for (int i = 0; i < N_L; i++) {
        std::cout << "\tL = " << L[i] << "... ";
        Clock2 clock(L[i], 2);
        clock.set_run_param(30000, 50000);
        auto binder = compute_binder(T, clock);
        data_clean2.insert_array(binder.data());
    }

    std::cout << "Performing 2D clock (2 spins) (disorder)\n";
    for (int i = 0; i < N_L; i++) {
        std::cout << "\tL = " << L[i] << "... ";
        Clock2 clock(L[i], 2);
        clock.set_run_param(30000, 50000);
        auto binder = compute_binder(T, clock, delta, N_run);
        data_disorder2.insert_array(binder.data());
    }

    std::cout << "Performing 3D clock (2 spins) (clean)\n";
    for (int i = 0; i < N_L; i++) {
        std::cout << "\tL = " << L[i] << "... ";
        Clock3 clock(L[i], 2);
        clock.set_run_param(30000, 50000);
        auto binder = compute_binder(T, clock);
        data_clean3.insert_array(binder.data());
    }

    std::cout << "Performing 3D clock (2 spins) (clean)\n";
    for (int i = 0; i < N_L; i++) {
        std::cout << "\tL = " << L[i] << "... ";
        Clock3 clock(L[i], 2);
        clock.set_run_param(30000, 50000);
        auto binder = compute_binder(T, clock, delta, N_run);
        data_disorder3.insert_array(binder.data());
    }

    std::ofstream of_clean2("binder_clean_clock2.dat"), of_disorder2("binder_disorder_clock2.dat"),
        of_clean3("binder_clean_clock3.dat"), of_disorder3("binder_disorder_clock3.dat");

    of_clean2 << data_clean2;
    of_clean3 << data_clean3;
    of_disorder2 << data_disorder2;
    of_disorder3 << data_disorder3;

    of_clean2.close();
    of_clean3.close();
    of_disorder2.close();
    of_disorder3.close();
}


/* test_xy()
 * Tests XY Model
 */
void test_xy(const std::array<double, N_pts> &T)
{
}
