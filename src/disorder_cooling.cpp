#include <iostream>
#include <array>
#include <random>
#include <type_traits>


/* compute_energy_clean()
 * Finds the energy of a clean model.
 */
template <typename TT, typename Model, size_t N>
std::array<TT, N> compute_energy_clean(const std::array<TT, N> &T, Model &model)
{
    if (!(std::is_same<double, TT>::value || std::is_same<float, TT>::value)) {
        std::cerr << "Error: Expected array of floar or double." << std::endl;
        exit(EXIT_FAILURE);
    } // Check for correct inputs.

    std::array<TT, N> E = {0};

    run_mc_energy(T, E, model);

    return E;
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
