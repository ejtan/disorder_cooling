#include <vector>
#include <iomanip>


template <int Row, int Col>
Data_matrix<Row, Col>::Data_matrix() : filled_row(0)
{
}


template <int Row, int Col>
Data_matrix<Row, Col>::Data_matrix(const Data_matrix &rhs) :
    filled_row(rhs.filled_row), data(rhs.data)
{
}


template <int Row, int Col>
void Data_matrix<Row, Col>::insert_data(const std::vector<double> &v)
{
    // Stores v[i] to data[i][filled_row]
    for (int i = 0; i < Row; i++)
        data[i * Col + filled_row] = v[i];

    // Incraments the number of rows filled
    filled_row++;
}


template <int Row, int Col>
std::ostream& operator<<(std::ostream &os, const Data_matrix<Row, Col> &rhs)
{
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Col; j++) {
            os << std::left << std::setw(12) << std::setprecision(8) << rhs.data[i * Col + j] << ' ';
        }
        os << '\n';
    } // Loop to write matrix
    os << std::flush;

    return os;
}

