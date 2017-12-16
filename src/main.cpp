#include <iostream>
#include <array>
#include <algorithm>

#include "../include/disorder_cooling.h"
#include "../include/ising2.h"
#include "../include/clock2.h"
#include "../include/xy2.h"


/*-------------------------------------------------------------------------------------------------
 * GLOBAL VARIABLES
 *-----------------------------------------------------------------------------------------------*/
const std::array<int, 3> L = {4, 8, 12};
const int N_pts = 60;
const double dT = 0.1;


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

    std::cout << "\nTesting Ising Model.\n";
    //test_ising(T);

    std::cout << "Testing Clock Model.\n";
    //test_clock(T);

    std::cout << "Testing XY Model.\n";
    test_xy(T);
}


/*-------------------------------------------------------------------------------------------------
 * FUNCTIONS
 *-----------------------------------------------------------------------------------------------*/
void test_ising(const std::array<double, N_pts> &T)
{
    double binder_dat[N_pts][3];

    for (int i = 0; i < 3; i++) {
        std::cout << "\tL = " << L[i] << "... ";
        Ising2 ising(L[i]);
        ising.set_run_param(30000, 50000);
        auto binder = compute_binder(T, ising);
        std::cout << "Done\n";

        for (int j = 0; j < N_pts; j++)
            binder_dat[j][i] = binder[j];
    }

    std::ofstream of("clean_ising_binder.dat");

    for (int i = 0; i < N_pts; i++) {
        of << T[i] << ' ';
        for (int j = 0; j < 3; j++)
            of << binder_dat[i][j] << ' ';
        of << '\n';
    }

    of.close();
}


void test_clock(const std::array<double, N_pts> &T)
{
    double binder_dat[N_pts][3];

    for (int i = 0; i < 3; i++) {
        std::cout << "\tL = " << L[i] << "... ";
        Clock2 clock(L[i], 2);
        clock.set_run_param(30000, 50000);
        auto binder = compute_binder(T, clock);
        std::cout << "Done\n";

        for (int j = 0; j < N_pts; j++)
            binder_dat[j][i] = binder[j];
    }

    std::ofstream of("clean_clock_binder.dat");

    for (int i = 0; i < N_pts; i++) {
        of << T[i] << ' ';
        for (int j = 0; j < 3; j++)
            of << binder_dat[i][j] << ' ';
        of << '\n';
    }

    of.close();
}


void test_xy(const std::array<double, N_pts> &T)
{
    double binder_dat[N_pts][3];

    for (int i = 0; i < 3; i++) {
        std::cout << "\tL = " << L[i] << "... ";
        XY2 xy(L[i]);
        xy.set_run_param(30000, 50000);
        auto binder = compute_binder(T, xy);
        std::cout << "Done\n";

        for (int j = 0; j < N_pts; j++)
            binder_dat[j][i] = binder[j];
    }

    std::ofstream of("clean_xy_binder.dat");

    for (int i = 0; i < N_pts; i++) {
        of << T[i] << ' ';
        for (int j = 0; j < 3; j++)
            of << binder_dat[i][j] << ' ';
        of << '\n';
    }

    of.close();
}
