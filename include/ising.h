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
        static const int warmup  = 30000;
        static const int measure = 500000;
        std::uniform_real_distribution<float> rand0;
        std::vector<int> spin;
        void sweep_lattice(const double beta, std::mt19937 &engine);

    public:
        Ising();
        Ising(const int L, const int dim);
        void init(const int L, const int dim);
        void set_spin();
        double sweep_energy(const double beta, std::mt19937 &engine);
};

#endif
