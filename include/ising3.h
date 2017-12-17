#ifndef ISING3_H
#define ISING3_H


#include <vector>
#include <random>

#include "model3.h"


/* class : Ising3
 * 3D Ising model class. Handles setting spin lattice and perfroming Monte Carlo sweeps using the
 * Metropolis Algorithm.
 */
class Ising3 : public Model3
{
    private:
        std::vector<int> spin;
        void sweep_lattice_clean(float beta, std::mt19937 &engine);
        void sweep_lattice_disorder(float beta, std::mt19937 &engine);

    public:
        Ising3() = default;
        Ising3(const int L);
        Ising3(const Ising3 &rhs);
        void set_spin();
        double sweep_energy(double beta, std::mt19937 &engine);
        double sweep_binder(double beta, std::mt19937 &engine);
};

#endif
