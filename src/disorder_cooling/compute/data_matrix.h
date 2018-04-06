#ifndef DATA_MATRIX_H
#define DATA_MATRIX_H

#include <iostream>
#include <vector>
#include <array>


// Forward Declarations
template <int Row, int Col> class Data_matrix;
template <int Row, int Col>
std::ostream& operator<<(std::ostream &os, const Data_matrix<Row, Col> &rhs);


/* Class: Data_matrix
 *
 * @template: Row = number of rows
 * @template: Col = number of columns
 *
 * Class which stores multiple columns of data in the format:
 * A_1,1   A_2,1   A_3,1   ...   A_N,1
 * A_1,2   A_2,2   A_3,2   ...   A_N,2
 * .       .       .       .     .
 * .       .       .        .    .
 * .       .       .         .   .
 * A_1,M   A_2,M   A_3_M   ...   A_N,M
 */
template <int Row, int Col>
class Data_matrix
{
    private:
        int filled_row;
        std::array<double, Row * Col> data;

    public:
        /* Constructors
         */
        Data_matrix();
        Data_matrix(const Data_matrix &rhs);

        /* insert_data()
         *
         * @input: v = data to insert
         *
         * Inserts the vector of data into the correct location.
         */
        void insert_data(const std::vector<double> &v);

        /* overloaded<<
         *
         * @input: os = ostream
         * @input: rhs = Data matrix to output.
         *
         * Outputs the matrix following the format outlined.
         */
        friend std::ostream& operator<< <Row, Col>(std::ostream &os, const Data_matrix &rhs);
};


#include "data_matrix.cpp"

#endif
