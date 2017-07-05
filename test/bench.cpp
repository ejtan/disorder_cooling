#include <iostream>
#include <array>
#include <algorithm>
#include <chrono>

#include "../include/disorder_cooling.h"


const double dT = 0.1;
const int n_pts = 200;
const int L     = 4;

void run_ising(const std::array<float, n_pts> &T)
{
    Ising ising(L, 2);

    auto start = std::chrono::system_clock::now();
    auto E = compute_energy_clean(T, ising);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "Time taken: " << elapsed_seconds.count() << "s." << std::endl;

    compute_entropy(E, T, "ising_clean.txt");
}


int main(int argc, char **argv)
{
    std::array<float, n_pts> T;
    int curr = 0;

    std::generate(T.begin(), T.end(), [&curr]() {
            curr++;
            return dT * static_cast<double>(curr);
    });

    run_ising(T);

    return 0;
}
