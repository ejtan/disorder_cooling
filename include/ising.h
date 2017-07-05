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
        std::vector<int> spin;
        void sweep_lattice(const double beta, std::mt19937 &engine);

    public:
        Ising(const int L, const int dim);
        Ising(const Ising &rhs);
        void init(const int L, const int dim);
        void set_spin();
        void set_exchange(const double delta);
        void set_exchange(const double J_val, const double p);
        double sweep_energy(const double beta, std::mt19937 &engine);
};


#endif
