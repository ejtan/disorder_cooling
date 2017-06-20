#ifndef CLOCK_H
#define CLOCK_H

#include <vector>
#include "model.h"


/* Clock Model class
 *
 * Class which handles Monte Carlo simulation of classical Clock Model.
 * Stores the value of cos into a table and is accessed using the spin array which
 * contains the indices to the cos table.
 */
class Clock : public Model
{
    private:
        int q;
        std::vector<int> spin;
        std::vector<double> cos_val;

        void sweep_lattice(const double beta, std::mt19937 &engine);
        double calculate_deltaE(const int pos, const int new_angle);
        double calculate_totalE();

    public:
        Clock(const int L, const int dim, const int _q);
        void init(const int L, const int dim, const int _q);
        void set_spin();
        void set_exchange(const double delta);
        void set_exchange(const double J_val, const double p);
        double sweep_energy(const double beta, std::mt19937 &engine);
};


#endif
