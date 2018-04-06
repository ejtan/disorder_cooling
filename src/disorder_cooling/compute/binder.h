#ifndef BINDER_H
#define BINDER_H

#include <vector>
#include <string>


/* compute_binder()
 *
 * @template: Model = model class
 *
 * @input: T = vector of temperatures used for computing correspnding binder ratios
 * @input: model = model we are simulating
 *
 * @output: vector<double> = energies corresponding to T
 *
 * Performs Monte Carlo simulation to compute binder ratios for a clean system.
 */
template <class Model>
std::vector<double> compute_binder(const std::vector<double> &T, Model &model);


/* compute_binder()
 *
 * @template: Model = model class
 *
 * @input: T = vector of temperatures used for computing correspnding binder ratios
 * @input: model = model we are simulating
 * @input: delta = width of disorder
 * @input: n_run = number of runs to perform
 *
 * @output: vector<double> = energies corresponding to T
 *
 * Performs Monte Carlo simulation to compute binder ratios for a clean system.
 */
template <class Model>
std::vector<double> compute_binder(const std::vector<double> &T, Model &model,
        double delta, int n_run);


#include "binder.cpp"

#endif
