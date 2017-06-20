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
        std::uniform_real_distribution<float> rand0;
        std::vector<int> spin;

        void sweep_lattice(const double beta, std::mt19937 &engine);

    public:
        XY();
        XY(const int _L, const int _dim);
        void init(const int _L, const int _dim);
        void set_spin();
        void set_exchange(const double delta);
        void set_exchange(const double J_val, const double p);
        double sweep_energy(const double beta, std::mt19937 &engine);
};

#endif
