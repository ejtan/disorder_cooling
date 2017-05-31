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
    private:
        int size;
        std::vector<int> spin;
        Neighbor_table neigh;
        Exchange_table J;

    public:
        Model();
};

#endif
