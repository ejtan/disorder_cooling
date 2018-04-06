#ifndef ENERGY_H
#define ENERGY_H

#include <vector>
#include <string>


/* compute_energy()
 *
 * @template: Model = model class
 *
 * @input: T = vector of temperatures used for computing corresponding energies
 * @input: model = model we are using
 *
 * @output: vector<double> = energies corresponding to temp
 *
 * Performs Monte Carlo simulation to calculate the energy for a clean system.
 */
template <class Model>
std::vector<double> compute_energy(const std::vector<double> &T, Model &model);


/* compute_energy()
 *
 * @template: Model = model class
 *
 * @input: T = vector of temperatures used for computing corresponding energies
 * @input: model = model we are using
 * @input: delta = width of disorder
 * @input: n_run = number of runs to perform
 *
 * @output: vector<double> = energies corresponding to temp
 *
 * Performs Monte Carlo simulation to calculate the energy for a disorder system.
 */
template <class Model>
std::vector<double> compute_energy(const std::vector<double> &T, Model &model,
        double delta, int n_run);


/* write_energy()
 *
 * @input: filename = name of file to output
 * @input: E = vector of energy
 * @input: T = vector of temperatures
 *
 * Writes E and T to a file. Format is: T E
 */
void write_energy(const std::string &filename,
        const std::vector<double> &E, const std::vector<double> &T);


/* write_energy()
 *
 * @input: filename = name of file to output
 * @input: E = vector of energy
 * @input: T = vector of temperatures
 * @input: N_spin = number of spins the model used.
 *
 * Computes the entropy and writes S and T to a file. Format is: T S.
 * Note that only N-1 points will be written due to integration.
 */
void write_entropy(const std::string &filename,
        const std::vector<double> &E, const std::vector<double> &T, int N_spin);


#include "energy.cpp"

#endif
