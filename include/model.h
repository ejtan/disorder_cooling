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
        int size, n_neigh;
        Neighbor_table neigh;
        Exchange_table J;

    public:
        Model();
        Model(const int L, const int dim);
        void init(const int L, const int dim);
        virtual void set_spin() = 0;
        virtual double sweep_energy(const double beta, std::mt19937 &engine) = 0;
};

#endif
