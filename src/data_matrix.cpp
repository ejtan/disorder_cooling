#include "../include/data_matrix.h"


/* Default constructor
 */
Data_matrix::Data_matrix() : N_row(0), N_col(0), filled_row(0), data(nullptr)
{
}


/* Constructor with arguments
 */
Data_matrix::Data_matrix(size_t r, size_t c) : N_row(r), N_col(c), filled_row(0)
{
    data = new double[r * c];
}


/* Copy constructor
 */
Data_matrix::Data_matrix(const Data_matrix &rhs) :
    N_row(rhs.N_row), N_col(rhs.N_col), filled_row(rhs.filled_row), data(rhs.data)
{
}
