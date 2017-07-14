#include <iostream>
#include <fstream>
#include <array>
#include <algorithm>
#include <random>
#include <type_traits>


/* compute_energy_clean()
 * Finds the energy of a clean model.
 */
template <typename TT, typename Model, size_t N>
std::array<TT, N> compute_energy(const std::array<TT, N> &T, Model &model)
{
    if (!(std::is_same<double, TT>::value || std::is_same<float, TT>::value)) {
        std::cerr << "Error: Expected array of floar or double." << std::endl;
        exit(EXIT_FAILURE);
    } // Check for correct inputs.

    std::array<TT, N> E = {0};

    run_mc_energy(T, E, model);

    return E;
}


/* compute_energy_disorder()
 * Finds the energy of a disorder model.
 */
template <typename TT, typename Model, size_t N>
std::array<TT, N> compute_energy(const std::array<TT, N> &T, Model &model,
        double delta, int n_run)
{
    if (!(std::is_same<double, TT>::value || std::is_same<float, TT>::value)) {
        std::cerr << "Error: Expected array of floar or double." << std::endl;
        exit(EXIT_FAILURE);
    } // Check for correct inputs.

    std::array<TT, N> E = {0};

    for (int run = 0; run < n_run; run++) {
        model.set_exchange(delta);
        run_mc_energy(T, E, model);
    } // Loop over runs

    // Normalize data
    std::transform(E.begin(), E.end(), E.begin(),
            [n_run](double val) { return val / static_cast<double>(n_run); });

    return E;
}


/* compute_entropy()
 * Takes in an energy and a tempearture array, computes the entropy, and outputs to a file.
 * There will be N-1 points in the output file due to the integration.
 */
template <typename TT, size_t N>
void compute_entropy(const std::array<TT, N> &E, const std::array<TT, N> &T,
        const std::string &filename)
{
    std::ofstream of(filename);
    const double ln2 = 0.693147180559945;

    std::array<TT, N> integral_values;

    for (size_t i = 0; i < N; i++)
        integral_values[i] = E[i] / (T[i] * T[i]);

    for (size_t i = 0; i < N - 1; i++)
        of << T[i] << ' ' << ln2 + (E[i] / T[i]) - trapezoid(T, integral_values, i) << '\n';

    of.close();
}


/*-------------------------------------------------------------------------------------------------
 * Intended Helper functions
 *-----------------------------------------------------------------------------------------------*/

/* run_mc_energy()
 * Performs Monte Carlo runs which computes the energy.
 */
template <typename TT, typename Model, size_t N>
void run_mc_energy(const std::array<TT, N> &T, std::array<TT, N> &E, Model model)
{
    std::random_device rd;
    std::mt19937 engine(rd());

    for (size_t i = 0; i < N; i++) {
        model.set_spin();
        E[i] += model.sweep_energy(1.0 / T[i], engine);
    } // Loop over all temperatures
}


/* trapezoid()
 * Perfroms integration with the trapezoidal rule with arbituary step sizes.
 */
template <typename TT, size_t N>
double trapezoid(const std::array<TT, N> &x, const std::array<TT, N> &y, int idx)
{
    double sum = 0.0;

    #pragma omp simd reduction(+:sum)
    for (size_t i = idx; i < N - 1; i++) {
        double dx = x[i + 1] - x[i];
        sum += dx * (y[i] + y[i + 1]);
    } // Loop to compute the integral

    return 0.5 * sum;
}
