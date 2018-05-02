#include <cmath>


/*-------------------------------------------------------------------------------------------------
 * PRIVATE METHODS
 *-----------------------------------------------------------------------------------------------*/
template <int Dim, int L, int q>
void Clock<Dim, L, q>::sweep_lattice_clean(float beta, std::mt19937 &engine)
{
    for (int i = 0; i < size; i++) {
        int pos = static_cast<int>(rand0(engine) * size);

        // Compute new angle
        int new_angle = static_cast<int>(rand0(engine) * q);
        while (new_angle == spin[pos])
            new_angle = static_cast<int>(rand0(engine) * q);

        float delta_E = 0.0;
        for (int i = 0; i < Dim * 2; i++) {
            int neigh_angle = spin[neigh[pos].neighbor[i]];
            int old_angle = spin[pos];

            int old_idx = (old_angle - neigh_angle + q) % q;
            int new_idx = (new_angle - neigh_angle + q) % q;

            delta_E += cos_val[old_idx] - cos_val[new_idx];
        } // Compute change in energy from all neighbors

        if (rand0(engine) < exp(-beta * delta_E))
            spin[pos] = new_angle;
    }
}


template <int Dim, int L, int q>
void Clock<Dim, L, q>::sweep_lattice_disorder(float beta, std::mt19937 &engine)
{
    for (int i = 0; i < size; i++) {
        int pos = static_cast<int>(rand0(engine) * size);

        // Compute new angle
        int new_angle = static_cast<int>(rand0(engine) * q);
        while (new_angle == spin[pos])
            new_angle = static_cast<int>(rand0(engine) * q);

        float delta_E = 0.0;
        for (int i = 0; i < Dim * 2; i++) {
            int neigh_angle = spin[neigh[pos].neighbor[i]];
            int old_angle = spin[pos];

            int old_idx = (old_angle - neigh_angle + q) % q;
            int new_idx = (new_angle - neigh_angle + q) % q;

            delta_E += J[pos].J_arr[i] * (cos_val[old_idx] - cos_val[new_idx]);
        } // Compute change in energy from all neighbors

        if (rand0(engine) < exp(-beta * delta_E))
            spin[pos] = new_angle;
    }
}


template <int Dim, int L, int q>
double Clock<Dim, L, q>::compute_energy()
{
    double E_tot = 0.0;

    for (int i = 0; i < size; i++) {
        int pos_angle = spin[i];
        int neigh1 = spin[neigh[i].neighbor[1]];
        int neigh2 = spin[neigh[i].neighbor[2]];

        int E_idx1 = (pos_angle - neigh1 + q) % q;
        int E_idx2 = (pos_angle - neigh2 + q) % q;

        if (isClean) {
            if (Dim == 2) {
                E_tot += -(cos_val[E_idx1] + cos_val[E_idx2]);
            } else if (Dim == 3) {
                int neigh3 = spin[neigh[i].neighbor[4]];
                int E_idx3 = (pos_angle - neigh3 + q) % q;

                E_tot += -(cos_val[E_idx1] + cos_val[E_idx2] + cos_val[E_idx3]);
            }
        } else {
            if (Dim == 2) {
                E_tot += -(J[i].J_arr[1] * cos_val[E_idx1] + J[i].J_arr[2] * cos_val[E_idx2]);
            } else if (Dim == 3) {
                int neigh3 = spin[neigh[i].neighbor[4]];
                int E_idx3 = (pos_angle - neigh3 + q) % q;

                E_tot += -(J[i].J_arr[1] * cos_val[E_idx1] + J[i].J_arr[2] * cos_val[E_idx2] +
                        J[i].J_arr[4] * cos_val[E_idx3]);
            }
        }
    }

    return E_tot;
}


template <int Dim, int L, int q>
double Clock<Dim, L, q>::compute_magnetization()
{
    double Mx = 0, My = 0;

    #pragma omp simd reduction(+:Mx, My)
    for (int i = 0; i < size; i++) {
        Mx += cos_val[spin[i]];
        My += sin_val[spin[i]];
    } // Loop over all sites

    return sqrt((Mx * Mx) + (My * My));
}

/*-------------------------------------------------------------------------------------------------
 * PUBLIC METHODS
 *-----------------------------------------------------------------------------------------------*/
template <int Dim, int L, int q>
Clock<Dim, L, q>::Clock()
{
    spin.resize(size);
    cos_val.resize(q);
    sin_val.resize(q);

    double dq = 2.0 * M_PI / static_cast<double>(q);

    std::generate(cos_val.begin(), cos_val.end(), [dq, i = 0]() mutable {
            double ret_val = cos(i * dq);
            i++;
            return ret_val;
    });
    std::generate(sin_val.begin(), sin_val.end(), [dq, i = 0]() mutable {
            double ret_val = cos(i * dq);
            i++;
            return ret_val;
    });
}


template <int Dim, int L, int q>
Clock<Dim, L, q>::Clock(const Clock &rhs) :
    Model<Dim, L>(rhs), spin(rhs.spin), cos_val(rhs.cos_val), sin_val(rhs.sin_val)
{
}


template <int Dim, int L, int q>
void Clock<Dim, L, q>::set_spin()
{
    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<int> dist(0, q);

    std::generate(spin.begin(), spin.end(), [&engine, &dist]()->int { return dist(engine); });
}


template <int Dim, int L, int q>
double Clock<Dim, L, q>::sweep_energy(double beta, std::mt19937 &engine)
{
    double E_tot = 0.0;

    if (isClean) {
        for (int i = 0; i < warmup; ++i)
            sweep_lattice_clean(beta, engine);

        for (int i = 0; i < measure; ++i) {
            sweep_lattice_clean(beta, engine);

            E_tot += compute_energy();
        } // measure sweeps
    } else {
        for (int i = 0; i < warmup; ++i)
            sweep_lattice_disorder(beta, engine);

        for (int i = 0; i < measure; ++i) {
            sweep_lattice_disorder(beta, engine);

            E_tot += compute_energy();
        } // measure sweeps
    } // Sweep clean or disrodered system

    return E_tot / static_cast<double>(measure * size);
}


template <int Dim, int L, int q>
double Clock<Dim, L, q>::sweep_binder(double beta, std::mt19937 &engine)
{
    double M2 = 0.0, M4 = 0.0;

    if (isClean) {
        for (int i = 0; i < warmup; ++i)
            sweep_lattice_clean(beta, engine);

        for (int i = 0; i < measure; ++i) {
            sweep_lattice_clean(beta, engine);

            double M = compute_magnetization();
            M2 += M * M;
            M4 += M * M * M * M;
        } // measure sweeps
    } else {
        for (int i = 0; i < warmup; ++i)
            sweep_lattice_disorder(beta, engine);

        for (int i = 0; i < measure; ++i) {
            sweep_lattice_disorder(beta, engine);

            double M = compute_magnetization();
            M2 += M * M;
            M4 += M * M * M * M;
        } // measure sweeps
    } // Sweep clean or disrodered system

    M2 /= static_cast<double>(measure);
    M4 /= static_cast<double>(measure);

    return 1.0 - (M4 / (3.0 * M2 * M2));
}
