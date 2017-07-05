#include <iostream>
#include <fstream>
#include <cstddef>
#include <cmath>
#include <array>
#include <algorithm>
#include <random>
#include <omp.h>


/* run_mc()
 *
 * Primary function for running the Monte Carlo simulation. Forks threads
 * to work on different parts of the temperature array.
 *
 * Note: for some reason, firstprivate cant copy if model is passed as reference.
 */
template<typename T, typename Model, size_t N>
void run_mc(const std::array<T, N> &temp, std::array<T, N> &E, Model model)
{
    int chunk;

    #pragma omp parallel shared(chunk) firstprivate(model)
    {
        #pragma omp single
        chunk =  N / omp_get_num_threads();

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
std::array<T, N> compute_energy_clean(const std::array<T, N> &temp, Model &model)
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
        Model &model, int n_run, double delta)
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
        Model &model, int n_run, double J, double P)
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


/* compute_entropy()
 *
 * Computes the entropy given arrays of energy E and temperature temp.
 * Uses traozeoidal rule on the entropy integral (set as a subroutine). Outputs data to a file.
 */
template <typename T, size_t N>
void compute_entropy(const std::array<T, N> &E, const std::array<T, N> &temp,
        const std::string &filename)
{
    std::ofstream of(filename);

    for (int i = 0; i < N - 1; i++)
        of << temp[i] << ' ' << log(2) + (E[i] / temp[i]) - trapezoid_entropy(E, temp, i) << '\n';

    of.close();
}


/* trapezoid_entropy()
 *
 * Performs integration with trapezoidal rule. Does not assume a constant dT, but computes
 * the step ever iteration.
 */
template <typename T, size_t N>
double trapezoid_entropy(const std::array<T, N> &E, const std::array<T, N> &temp, int idx)
{
    double sum = 0.0;

    for (int i = idx; i < N - 1; i++) {
        double dT = temp[i + 1] - temp[i];
        sum += dT * ((E[i] / (temp[i] * temp[i])) + (E[i + 1] / (temp[i + 1] * temp[i + 1])));
    } // Loop to compute the integral

    return sum * 0.5;
}
