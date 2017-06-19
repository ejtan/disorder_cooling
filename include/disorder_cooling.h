#ifndef DISORDER_COOLING_H
#define DISORDER_COOLING_H

#include "../include/ising.h"
#include <array>


/* Header file for Monte Carlo simulations of Classica spin models.
 *
 * Contains methods to generate data and output to files.
 * Implmentation in header file because C++ is fucking dumb.
 */

/* run_mc()
 *
 * Primary function for running the Monte Carlo simulation. Forks threads
 * to work on different parts of the temperature array.
 */
template<typename T, typename Model, size_t N>
void run_mc(const std::array<T, N> &temp, std::array<T, N> &E, Model model)
{
    int chunk;

    // TODO:
    // At some point, figure out how to split uneven work among threads
    #pragma omp paralle shared(chunk, temp) threadprivate(model) num_threads(4)
    {
        #pragma omp single
        {
            chunk =  N / omp_get_num_threads();
        } // Single Region

        int thd_id = omp_get_thread_num();
        std::random_device rd;
        std::mt19937 engine(rd());

        for (int i = chunk * thd_id; i < (chunk * (thd_id + 1)); i++) {
            model.set_spin();
            E[i] += model.sweep_energy(1.0 / temp[i], engine);
        } // Loop over chunk of work
    } // Parallel Region
}


/* compute_energy_clean()
 *
 * Computes the energy of a given clean model.
 */
template<typename T, typename Model, size_t N>
std::array<T, N> compute_energy_clean(const std::array<T, N> &temp, Model model)
{
    if (!(std::is_same<double, T>::value || std::is_same<float, T>::value)) {
        std::cerr << "Error: Expected array of float or double" << std::endl;
        exit(EXIT_FAILURE);
    } // Check Model Inputs

    std::array<T, N> E;

    E.fill(0);
    run_mc(temp, E, model);

    return E;
}


/* compute_energy_disorder()
 *
 * Computes the energy of a given model with bond disorder defined in the
 * initalzation of the model. Method does continuous distribution.
 * Returns array with corresponding energies.
 */
template<typename T, typename Model, size_t N>
std::array<T, N> compute_energy_disorder(const std::array<T, N> &temp,
        Model model, const int n_run, const double delta)
{
    if (!(std::is_same<double, T>::value || std::is_same<float, T>::value)) {
        std::cerr << "Error: Expected array of float or double" << std::endl;
        exit(EXIT_FAILURE);
    } // Check Input

    std::array<double, N> E;

    // Fill Array to 0
    E.fill(0);

    for (int run = 0; run < n_run; run++) {
        model.set_exchange(delta);
        run_mc(temp, E, model);
    } // Loop over runs

    // Normalize data
    std::transform(E.begin(), E.end(), E.begin(),
            [n_run](double val) { return val / static_cast<double>(n_run);
    });

    return E;
}


/* compute_energy_disorder()
 *
 * Overloaded function for computing the energy of a given model with
 * discrete bond disorder
 */
template<typename T, typename Model, size_t N>
std::array<T, N> compute_energy_disorder(const std::array<T, N> &temp,
        Model model, const int n_run, const double J, const double P)
{
    if (!(std::is_same<double, T>::value || std::is_same<float, T>::value)) {
        std::cerr << "Error: Expected array of float or double" << std::endl;
        exit(EXIT_FAILURE);
    } // Check Model Inputs

    std::array<T, N> E;

    E.fill(0);

    for (int run = 0; run < n_run; run++) {
        model.set_exchange(J, P);
        run_mc(temp, E, model);
    }

    // Normalize data
    std::transform(E.begin(), E.end(), E.begin(),
            [n_run](double val) { return val / static_cast<double>(n_run);
    });

    return E;
}


#endif
