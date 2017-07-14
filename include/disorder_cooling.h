#ifndef DISORDER_COOLING_H
#define DISORDER_COOLING_H


#include <array>
#include <string>

#include "ising2.h"


/* Header file for Monte Carlo simulations of classical spin models.
 *
 * Provides functions for running the simulation, gathering data, and processing data to a file.
 */
template <typename TT, typename Model, size_t N>
std::array<TT, N> compute_energy_clean(const std::array<TT, N> &T, Model &model);

template <typename TT, size_t N>
void compute_entropy(const std::array<TT, N> &E, const std::array<TT, N> &T,
        const std::string &filename);


/*-------------------------------------------------------------------------------------------------
 * Intended Helper functions
 *-----------------------------------------------------------------------------------------------*/

template <typename TT, typename Model, size_t N>
void run_mc_energy(const std::array<TT, N> &T, std::array<TT, N> &E, Model model);

template <typename TT, size_t N>
double trapezoid(const std::array<TT, N> &x, const std::array<TT, N> &y, int idx);


#include "../src/disorder_cooling.cpp"

#endif
