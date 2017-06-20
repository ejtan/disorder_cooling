#ifndef XY_H
#define XY_H

#include "clock.h"

/* XY Model Class
 *
 * Class which handles Monte Carlo simulations for the classical XY Model.
 * Uses the clock model as a base and sets the number of discrete spins to 50.
 * At this number of spins, the model converges from the clock to XY model.
 */
class XY : public Clock
{
    private:
        std::vector<int> spin;

    public:
        XY();
        XY(const int L, const int dim);
        void init(const int L, const int dim);
        void set_spin();
        void set_exchange(const double delta);
        void set_exchange(const double J_val, const double p);
        double sweep_energy(const double beta, std::mt19937 &engine);
};

#endif
