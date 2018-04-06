#ifndef UTILITY_H
#define UTILITY_H


/* run_mc_energy()
 *
 * @template: Model = model class
 *
 * @input: T = vector of temperatures
 * @input: E = output vector of energies
 * @input: model = Model we are using
 *
 * Performs Monte Carlo simulations to compute E[i] corresponding to T[i].
 */
template <class Model>
void run_mc_energy(const std::vector<double> &T, std::vector<double> &E, Model model);


/* run_mc_energy()
 *
 * @template: Model = model class
 *
 * @input: T = vector of temperatures
 * @input: binder = output vector of binder ratios
 * @input: model = Model we are using
 *
 * Performs Monte Carlo simulations to compute binder[i] corresponding to T[i].
 */
template <class Model>
void run_mc_binder(const std::vector<double> &T, std::vector<double> &binder, Model model);


/* trapezoid()
 *
 * @input: x = values on the x-axis
 * @input: y = f(x[i])
 * @input: idx = starting index of integration
 *
 * Performs integration from x[idx] to x[N-1] using trapezoidal method.
 */
double trapezoid(const std::vector<double> &x, const std::vector<double> &y, int idx);


#include "utility.cpp"

#endif
