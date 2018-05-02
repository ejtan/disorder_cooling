#ifndef CLOCK_H
#define CLOCK_H

#include <vector>
#include <random>
#include "model.h"


/* Class: Clock
 *
 * @template: Dim = number of Dimensions
 * @template: L = number of sites per dimension
 * @template: q = number of possiable spins
 *
 * Clock Model class which handles setting the lattice and performing Monte Carlo sweeps with the
 * Metropolis Algorithm.
 */
template <int Dim, int L, int q>
class Clock : public Model<Dim, L>
{
    using Model<Dim, L>::size;
    using Model<Dim, L>::measure;
    using Model<Dim, L>::warmup;
    using Model<Dim, L>::neigh;
    using Model<Dim, L>::rand0;
    using Model<Dim, L>::J;
    using Model<Dim, L>::isClean;

    private: // Variables
        std::vector<int> spin;
        std::vector<double> cos_val, sin_val;

    private: // Methods
        /* sweep_lattice_clean()
         *
         * @param: beta = 1 / T
         * @param: engine = random number generator engine
         *
         * Perfroms clean lattice sweep. Goes through L ^ Dim elements.
         */
        void sweep_lattice_clean(float beta, std::mt19937 &engine);

        /* sweep_lattice_disorder()
         *
         * @param: beta = 1 / T
         * @param: engine = random number generator engine
         *
         * Perfroms disordered lattice sweeps. Goes through L ^ Dim elements.
         */
        void sweep_lattice_disorder(float beta, std::mt19937 &engine);

        /* compute_energy()
         *
         * Computes the energy of the lattice.
         */
        double compute_energy();

        /* compute_magnetization()
         *
         * Computes the magnetization of the lattice.
         */
        double compute_magnetization();

    public:
        /* Constructors
         */
        Clock();
        Clock(const Clock &rhs);

        /* Implamentation of virtual functions
         * See modeh.h for description
         */
        void set_spin();
        double sweep_energy(double beta, std::mt19937 &engine);
        double sweep_binder(double beta, std::mt19937 &engine);
};

#include "clock.cpp"

#endif
