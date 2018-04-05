#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <random>
#include "../tables/neighbor.h"
#include "../tables/exchange.h"


/* class: Model
 *
 * @template: Dim = number of dimensions
 * @template: L = number of lattice sites in a dimension
 *
 * Base class for 2D classical spin models.
 */
template <int Dim, int L>
class Model
{
    static_assert(Dim == 2 || Dim == 3, "Error: Dim must be 2 or 3.");
    static_assert(L > 0, "Error: L must be positive and nonzero.");

    protected:
        int warmup = 30000;
        int measure = 500000;
        int size;
        bool isClean;
        std::uniform_real_distribution<float> rand0;
        std::vector< Neighbor<Dim, L> > neigh;
        std::vector< Exchange<Dim> > J;

    public:
        /* Constructors
         */
        Model();
        Model(const Model &rhs);

        /* virtual set_spin()
         *
         * Initalizes the spins in the system.
         */
        virtual void set_spin() = 0;

        /* virtual sweep_energy()
         *
         * @param: beta = 1 / T
         * @param: engine = mt19937 random number generator
         *
         * @return: double representing the energy at beta.
         */
        virtual double sweep_energy(double beta, std::mt19937 &engine) = 0;

        /* virtual sweep_binder()
         *
         * @param: beta = 1 / T
         * @param: engine = mt19937 random number generator
         *
         * @return: double representing the binder ratio at beta.
         */
        virtual double sweep_binder(double beta, std::mt19937 &engine) = 0;

        /* set_exchange()
         *
         * @param: delta = width of the uniform distribution to pick disorders from, centered at 1.
         *
         * Sets the bond disorder.
         */
        void set_exchange(double delta);

        /* set_run_param()
         *
         * @param: w = warmup steps
         * @param: m = measure steps
         *
         * Overwrites the default initalized parameters.
         */
        void set_run_param(int w, int m);
};


#include "model.cpp"

#endif
