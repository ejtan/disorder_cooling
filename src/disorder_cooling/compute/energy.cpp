#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <cmath>

#include "utility.h"


template <class Model>
std::vector<double> compute_energy(const std::vector<double> &T, Model &model)
{
    std::vector<double> E(T.size());
    std::fill(E.begin(), E.end(), 0);
    run_mc_energy(T, E, model);

    return E;
}


template <class Model>
std::vector<double> compute_energy(const std::vector<double> &T, Model &model,
        double delta, int n_run)
{
    std::vector<double> E(T.size());
    std::fill(E.begin(), E.end(), 0);

    for (int run = 0; run < n_run; run++) {
        model.set_exchange(delta);
        run_mc_energy(T, E, model);
    } // Loop over runs

    // Normalize data
    std::transform(E.begin(), E.end(), E.begin(),
            [n_run](double val)->double { return val / static_cast<double>(n_run); });

    return E;
}


void write_energy(const std::string &filename,
        const std::vector<double> &E, const std::vector<double> &T)
{
    std::ofstream of(filename);
    int N = T.size();

    of << "#Temperature Energy\n";
    for (int i = 0; i < N; i++)
        of << std::left << std::setw(12) << std::setprecision(8) << T[i] << ' ' << E[i] << '\n';
    of << std::flush;

    of.close();
}


void write_entropy(const std::string &filename,
        const std::vector<double> &E, const std::vector<double> &T, int N_spin)
{
    std::ofstream of(filename);
    std::vector<double> integration_vals(T.size());
    int N = T.size();

    for (int i = 0; i < N; i++)
        integration_vals[i] = E[i] / (T[i] * T[i]);

    of << "#Temperature Entropy\n";
    for (int i = 0; i < N - 1; i++)
        of << std::setw(12) << std::setprecision(8) << T[i] << " "
            << log(N_spin) + (E[i] / T[i]) - trapezoid(T, integration_vals, i) << '\n';
}
