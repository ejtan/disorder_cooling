#ifndef ISING_H
#define ISING_H

#include <vector>
#include "model.h"


/* Ising Model class
 *
 * Class which handles Monte Carlo simulation of classical Ising Model.
 */
class Ising : public Model
{
    private:
        std::uniform_real_distribution<float> rand0;
        std::vector<int> spin;
        void sweep_lattice(const double beta, std::mt19937 &engine);

    public:
        Ising();
        Ising(const int _L, const int _dim);
        void init(const int _L, const int _dim);
        void set_spin();
        void set_exchange(const double delta);
        void set_exchange(const double J_val, const double p);
        double sweep_energy(const double beta, std::mt19937 &engine);
};

#endif
