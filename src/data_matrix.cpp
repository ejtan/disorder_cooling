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


/* Destuctor
 */
Data_matrix::~Data_matrix()
{
    if (data)
        delete[] data;
}


/* insert_array()
 * Takes a pointer to an array of data and inserts data. Assumes that the
 * input has N_row amount of data. Note that the data is stored in column major ordering.
 */
void Data_matrix::insert_array(const double *input)
{
    int offset = filled_row * N_row;

    for (size_t i = 0; i < N_row; i++)
        data[i * N_col + offset] = input[i];

    filled_row++;
}
