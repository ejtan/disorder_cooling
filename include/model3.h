#ifndef MODEL3_H
#define MODEL3_H


#include <vector>
#include <random>

#include "neighbor.h"
#include "exchange.h"


/* Base class for 3D Classical spin models.
 */
class Model3
{
    protected:
        size_t warmup  = 30000;
        size_t measure = 500000;
        static const int n_neigh = 6;
        size_t size;
        bool isClean;
        std::uniform_real_distribution<float> rand0;
        std::vector<Neighbor<3>> neigh;
        std::vector<Exchange<3>> J;

    public:
        Model3();
        Model3(const int L);
        Model3(const Model3 &rhs);
        virtual void set_spin() = 0;
        virtual double sweep_energy(double beta, std::mt19937 &engine) = 0;
        virtual double sweep_binder(double beta, std::mt19937 &engine) = 0;
        void set_exchange(double delta);
        std::vector<Exchange<3>> get_exchange() const;
        void set_run_param(size_t Warmup, size_t Measure);
};

#endif
