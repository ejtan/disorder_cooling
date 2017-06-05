#include "../include/model.h"


/* Default Constructor
 */
Model::Model()
{
}


/* Constructor with arguments.
 */
Model::Model(const int L, const int dim)
{
    size = L * L;
    neigh.init(L, dim);
    J.init(L, dim);
    J.generate_clean();
}


/* init()
 *
 * Behaves just as the constructor with arguments.
 */
void Model::init(const int L, const int dim)
{
    size = L * L;
    neigh.init(L, dim);
    J.init(L, dim);
    J.generate_clean();
}
