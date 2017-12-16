#ifndef CLOCK2_H
#define CLOCK2_H


#include <vector>
#include <random>
#include "model2.h"


/* Class : Clock2
 * Class for handeling Monte Carlo simulation of a clock model. Stores the values of the spin
 * (we only need the cosine values) into a table to be accessed by an index in the spin array.
 */
class Clock2 : public Model2
{
    private:
        int q;
        std::vector<int> spin;
        std::vector<double> cos_val;
        std::vector<double> sin_val;

        void sweep_lattice_clean(float beta, std::mt19937 &engine);
        void sweep_lattice_disorder(float beta, std::mt19937 &engine);

    public:
        Clock2() = default;
        Clock2(const int L, const int _q);
        Clock2(const Clock2 &rhs);
        void set_spin();
        double sweep_energy(double beta, std::mt19937 &engine);
        double sweep_binder(double beta, std::mt19937 &engine);
};

#endif
