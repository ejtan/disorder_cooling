#include <cmath>


/*-------------------------------------------------------------------------------------------------
 * PRIVATE METHODS
 *-----------------------------------------------------------------------------------------------*/
template <int Dim, int L>
void Ising<Dim, L>::sweep_lattice_clean(float beta, std::mt19937 &engine)
{
    for (int i = 0; i < size; i++) {
        int pos = static_cast<int>(rand0(engine) * size);
        float neigh_sum = 0.0;

        for (int j = 0; j < 2 * Dim; j++)
            neigh_sum += spin[neigh[pos].neighbor[j]];

        float delta_E = 2.0 * spin[pos] * neigh_sum;

        // Accept / reject flip
        if (rand0(engine) < exp(-beta * delta_E))
            spin[pos] = -spin[pos];
    } // Sweep over sites
}


template <int Dim, int L>
void Ising<Dim, L>::sweep_lattice_disorder(float beta, std::mt19937 &engine)
{
    for (int i = 0; i < size; i++) {
        int pos = static_cast<int>(rand0(engine) * size);
        float neigh_sum = 0.0;

        for (int j = 0; j  < 2 * Dim; j++)
            neigh_sum += J[pos].J_arr[j] * spin[neigh[pos].neighbor[j]];

        float delta_E = 2.0 * spin[pos] * neigh_sum;

        // Accept / reject flip
        if (rand0(engine) < exp(-beta * delta_E))
            spin[pos] = -spin[pos];
    } // Sweep over sites
}


/*-------------------------------------------------------------------------------------------------
 * PUBLIC METHODS
 *-----------------------------------------------------------------------------------------------*/
template <int Dim, int L>
Ising<Dim, L>::Ising()
{
    spin.resize(size);
}


template <int Dim, int L>
Ising<Dim, L>::Ising(const Ising &rhs) : Model<Dim, L>(rhs), spin(rhs.spin)
{
}


template <int Dim, int L>
void Ising<Dim, L>::set_spin()
{
    for (auto &&s : spin)
        s = 1;
}


template <int Dim, int L>
double Ising<Dim, L>::sweep_energy(double beta, std::mt19937 &engine)
{
    double E_tot = 0.0;

    if (isClean) {
        for (int i = 0; i < warmup; i++)
            sweep_lattice_clean(beta, engine);

        for (int i = 0; i < measure; i++) {
            sweep_lattice_clean(beta, engine);

            if (Dim == 2) {
                #pragma omp simd reduction(+:E_tot)
                for (int j = 0; j < size; j++)
                    E_tot += -spin[j] * (spin[neigh[j].neighbor[0]] + spin[neigh[j].neighbor[1]]);
            } else if (Dim == 3) {
                #pragma omp simd reduction(+:E_tot)
                for (int j = 0; j < size; j++)
                    E_tot += -spin[j] * (spin[neigh[j].neighbor[0]] + spin[neigh[j].neighbor[1]] +
                                         spin[neigh[j].neighbor[4]]);
            } // Measure energy
        } // Measure sweeps
    } else {
        for (int i = 0; i < warmup; i++)
            sweep_lattice_disorder(beta, engine);

        for (int i = 0; i < measure; i++) {
            sweep_lattice_disorder(beta, engine);

            if (Dim == 2)
                #pragma omp simd reduction(+:E_tot)
                for (int j = 0; j < size; j++)
                    E_tot += -spin[j] * (J[j].J_arr[0] * spin[neigh[j].neighbor[0]] +
                                         J[j].J_arr[1] * spin[neigh[j].neighbor[1]]);
            else if (Dim == 3)
                #pragma omp simd reduction(+:E_tot)
                for (int j = 0; j < size; j++)
                    E_tot += -spin[j] * (J[j].J_arr[0] * spin[neigh[j].neighbor[0]] +
                                         J[j].J_arr[1] * spin[neigh[j].neighbor[1]] +
                                         J[j].J_arr[3] * spin[neigh[j].neighbor[4]]);
        } // Measure sweeps
    } // Sweep clean or disorder system

    return E_tot / static_cast<double>(measure * size);
}


template <int Dim, int L>
double Ising<Dim, L>::sweep_binder(double beta, std::mt19937 &engine)
{
    double M2 = 0, M4 = 0;

    if (isClean) {
        for (int i = 0; i < warmup; i++)
            sweep_lattice_clean(beta, engine);

        for (int i = 0; i < measure; i++) {
            sweep_lattice_clean(beta, engine);

            double M = 0.0;
            #pragma omp simd reduction(+:M)
            for (int j = 0; j < size; j++)
                M += spin[j];

            M2 += M * M;
            M4 += M * M * M * M;
        } // Measure sweeps
    } else {
        for (int i = 0; i < warmup; i++)
            sweep_lattice_disorder(beta, engine);

        for (int i = 0; i < measure; i++) {
            sweep_lattice_disorder(beta, engine);

            double M = 0.0;
            #pragma omp simd reduction(+:M)
            for (int j = 0; j < size; j++)
                M += spin[j];

            M2 += M * M;
            M4 += M * M * M * M;
        } // Measure sweeps
    } // Sweep clean or disorder system

    // Normalize data
    M2 /= static_cast<double>(measure);
    M4 /= static_cast<double>(measure);

    return 1.0 - (M4 / (3.0 * M2 * M2));
}
