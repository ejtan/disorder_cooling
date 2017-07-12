#ifndef MODEL2_H
#define MODEL2_H


#include <vector>
#include <random>
#include "neighbor.h"


/* Base class for 2D Classical spin models.
 */
class Model2
{
    protected:
        static const int warmup  = 30000;
        static const int measure = 500000;
        int size;
        std::uniform_real_distribution<float> rand0;
        std::vector<Neighbor<2>> neigh;
        std::vector<double> J;
    public:
        Model2();
        Model2(const int L);
        Model2(const Model2 &rhs);
};

#endif
