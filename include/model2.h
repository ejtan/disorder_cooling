#ifndef MODEL2_H
#define MODEL2_H


#include <vector>
#include <random>
#include "neighbor.h"
#include "exchange.h"


/* Base class for 2D Classical spin models.
 */
class Model2
{
    protected:
        static const int warmup  = 30000;
        static const int measure = 500000;
        static const int n_neigh = 4;
        int size;
        bool isClean;
        std::uniform_real_distribution<float> rand0;
        std::vector<Neighbor<2>> neigh;
        std::vector<Exchange<2>> J;

    public:
        Model2();
        Model2(const int L);
        Model2(const Model2 &rhs);
        virtual void set_spin() = 0;
        virtual double sweep_energy(double beta, std::mt19937 &engine) = 0;
        void set_exchange(double delta);
        std::vector<Exchange<2>> get_exchange() const;
};

#endif
