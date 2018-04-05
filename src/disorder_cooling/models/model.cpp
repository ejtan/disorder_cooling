template <int Dim, int L>
Model<Dim, L>::Model() : isClean(true), rand0(0.0, 1.0)
{
    if (Dim == 2)
        size = L * L;
    else if (Dim == 3)
        size = L * L * L;

    neigh.resize(size);
    J.resize(size);

    for (int i = 0; i < size; ++i)
        neigh[i].set_neighbors(i);
}


template <int Dim, int L>
Model<Dim, L>::Model(const Model &rhs) :
    size(rhs.size), isClean(rhs.isClean), neigh(rhs.neigh), J(rhs.J)
{
}



template <int Dim, int L>
void Model<Dim, L>::set_exchange(double delta)
{
    if (isClean)
        isClean = false;

    std::random_device rd;
    std::mt19937 engine(rd());
    float r_val, J_val;

    for (int i = 0; i < size; i++) {
        // 0 - 2 bond
        r_val = rand0(engine);
        if (rand0(engine) > 0.5) J_val = 1.0 - (delta * r_val / 2.0);
        else                     J_val = 1.0 + (delta * r_val / 2.0);
        J[i].J_arr[0]                    = J_val;
        J[neigh[i].neighbor[0]].J_arr[2] = J_val;

        // 1 - 3 bond
        r_val = rand0(engine);
        if (rand0(engine) > 0.5) J_val = 1.0 - (delta * r_val / 2.0);
        else                     J_val = 1.0 + (delta * r_val / 2.0);
        J[i].J_arr[1]                    = J_val;
        J[neigh[i].neighbor[1]].J_arr[3] = J_val;

        if (Dim == 3) {
            // 4 - 5 bond
            r_val = rand0(engine);
            if (rand0(engine) > 0.5) J_val = 1.0 - (delta * r_val / 2.0);
            else                     J_val = 1.0 + (delta * r_val / 2.0);
            J[i].J_arr[4]                    = J_val;
            J[neigh[i].neighbor[4]].J_arr[5] = J_val;
        } // Set 3D bond if needed
    } // Set bond
}


template <int Dim, int L>
void Model<Dim, L>::set_run_param(int w, int m)
{
    warmup = w;
    measure = m;
}
