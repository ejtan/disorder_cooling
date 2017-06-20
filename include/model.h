#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <random>
#include "neighbor_table.h"
#include "exchange_table.h"


/* Base class for Classical spin models.
 */
class Model
{
    protected:
        static const int warmup  = 30000;
        static const int measure = 500000;
        int size, n_neigh;
        std::uniform_real_distribution<float> rand0;
        Neighbor_table neigh;
        Exchange_table J;

    public:
        Model();
        Model(const int L, const int dim);
        void init(const int L, const int dim);
        virtual void set_spin() = 0;
        void set_exchange(const double delta);
        void set_exchange(const double J_val, const double p);
        virtual double sweep_energy(const double beta, std::mt19937 &engine) = 0;
};


#endif
