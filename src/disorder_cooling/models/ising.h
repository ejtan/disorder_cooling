#ifndef ISING_H
#define ISING_H

#include <vector>
#include <random>
#include "model.h"


/* Class: Ising
 *
 * @template: Dim = number of dimensions
 * @template: L = number of sites per dimension
 *
 * Ising class which handles setting the lattice and perfroming Monte Carlo sweeps using the
 * Metropolis Algorithm.
 */
template <int Dim, int L>
class Ising : public Model<Dim, L>
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

    private: // Functions
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

    public:
        /* Constructors
         */
        Ising();
        Ising(const Ising &rhs);

        /* Implamentation of virtual functions
         * See model.h for description
         */
        void set_spin();
        double sweep_energy(double beta, std::mt19937 &engine);
        double sweep_binder(double beta, std::mt19937 &engine);
};


#include "ising.cpp"

#endif
