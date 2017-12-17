#ifndef DATA_MATRIX_H
#define DATA_MATRIX_H


/* Data_matrix class
 * Matrix for storing multiple columns of data.
 * Provides only methods for inserting data and outputing with an ostream operator.
 */
class Data_matrix
{
    private:
        size_t Nrow, Ncol;
        int filled_row;
        double *data;

    public:
        Data_matrix();
        Data_matrix(size_t r, size_t c);
        Data_matrix(const Data_matrix &rhs);
        ~Data_matrix();
};

#endif
