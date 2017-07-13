#ifndef ISING2_H
#define ISING2_H


#include <vector>
#include <random>
#include "model2.h"


/* class : Ising2
 * 2D Ising model class. Handles setting spin lattice and perfroming Monte Carlo sweeps using the
 * Metropolis Algorithm.
 */
class Ising2 : public Model2
{
    private:
        std::vector<int> spin;
        void sweep_lattice_clean(float beta, std::mt19937 &engine);

    public:
        Ising2() = default;
        Ising2(const int L);
        Ising2(const Ising2 &rhs);
        void set_spin();
        double sweep_energy(double beta, std::mt19937 &engine);

};

#endif
