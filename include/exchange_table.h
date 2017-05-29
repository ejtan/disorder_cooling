#ifndef EXCHANGE_TABLE_H
#define EXCHANGE_TABLE_H

#include <vector>
#include "neighbor_table.h"


/* Class : Exchange_table
 *
 * Table which set the exchange bond between neighboring sites.
 * Utilizes a neighbor table to set the values. Can take both discrete and
 * continuous distrubitons.
 *
 * If discrete, than +/- J values are set with probability p for +j and (1-p)
 * for -J.
 *
 * A continuous distribution set a box centered at J = 1 and ranges from
 * J = [1 + delta / 2, 1 + delta / 2].
 */
class Exchange_table
{
    private:
        int size, n_neigh;
        Neighbor_table neigh;
        std::vector<double> table;

    public:
        Exchange_table();
        Exchange_table(const int L, const int dim);
        void init(const int L, const int dim);
        void generate_discrete(const double J, const double prob);
        void generate_continuous(const double delta);
        double operator[](const int idx) const;
        std::vector<double> get_table();
};

#endif
