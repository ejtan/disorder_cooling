#ifndef CLOCK2_H
#define CLOCK2_H


#include <vector>
#include <random>

#include "model3.h"


/* Class : Clock3
 * Class for handeling Monte Carlo simulation of a clock model. Stores the values of the spin
 * (we only need the cosine values) into a table to be accessed by an index in the spin array.
 */
class Clock3 : public Model3
{
    private:
        int q;
        std::vector<int> spin;
        std::vector<double> cos_val;

        void sweep_lattice_clean(float beta, std::mt19937 &engine);
        void sweep_lattice_disorder(float beta, std::mt19937 &engine);

    public:
        Clock3() = default;
        Clock3(const int L, const int _q);
        Clock3(const Clock3 &rhs);
        void set_spin();
        double sweep_energy(double beta, std::mt19937 &engine);
};

#endif
