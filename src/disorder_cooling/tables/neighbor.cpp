#include <type_traits>


template <int Dim, int L>
void Neighbor<Dim, L>::set_neighbors(int pos)
{
    static_assert(Dim == 2 || Dim == 3, "Dim must be 2 or 3.");

    int N2 = L * L; // Number of sites on a 2D layer

    // 0-th neighbor
    if ((pos % N2) - L < 0)   neighbor[0] = pos + N2 - L;
    else                      neighbor[0] = pos - L;
    // 1-st neighbor
    if ((pos + 1) % L == 0)   neighbor[1] = pos + 1 - L;
    else                      neighbor[1] = pos + 1;
    // 2-nd neighbor
    if ((pos + L) % N2 < L)   neighbor[2] = pos + L - N2;
    else                      neighbor[2] = pos + L;
    // 3-rd neighbor
    if (pos % L == 0)         neighbor[3] = pos + L - 1;
    else                      neighbor[3] = pos - 1;

    if (Dim == 3) {
        int N3 = L * L * L; // Total number of sites in a 3D lattice

        // 4-th neighbor
        if (pos - N2 < 0)   neighbor[4] = pos + N3 - N2;
        else                neighbor[4] = pos - N2;
        // 5-th neighbor
        if (pos + N2 >= N3) neighbor[5] = pos - N3 + N2;
        else                neighbor[5] = pos + N2;
    } // Set 3D layer (if needed)
}
