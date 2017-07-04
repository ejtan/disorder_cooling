#ifndef DISORDER_COOLING_H
#define DISORDER_COOLING_H

#include "../include/ising.h"
#include "../include/clock.h"
#include "../include/xy.h"

#include <array>
#include <cstddef>


/* Header file for Monte Carlo simulations of Classica spin models.
 *
 * Contains methods to generate data and output to files.
 * Implementation in .cpp file included at the bottom due to dumb reasons.
 */

template<typename T, typename Model, size_t N>
void run_mc(const std::array<T, N> &temp, std::array<T, N> &E, Model model);

template<typename T, typename Model, size_t N>
std::array<T, N> compute_energy_clean(const std::array<T, N> &temp, Model model);

template<typename T, typename Model, size_t N>
std::array<T, N> compute_energy_disorder(const std::array<T, N> &temp,
        Model model, const int n_run, const double delta);

template<typename T, typename Model, size_t N>
std::array<T, N> compute_energy_disorder(const std::array<T, N> &temp,
        Model model, const int n_run, const double J, const double P);

#include "../src/disorder_cooling.cpp"

#endif
