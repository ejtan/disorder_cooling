template <int Dim, int L>
XY<Dim, L>::XY(const XY &rhs) : Clock<Dim, L, q>(rhs)
{
}


template <int Dim, int L>
void XY<Dim, L>::set_spin()
{
    Clock<Dim, L, q>::set_spin();
}


template <int Dim, int L>
double XY<Dim, L>::sweep_energy(double beta, std::mt19937 &engine)
{
    return Clock<Dim, L, q>::sweep_energy(beta, engine);
}


template <int Dim, int L>
double XY<Dim, L>::sweep_binder(double beta, std::mt19937 &engine)
{
    return Clock<Dim, L, q>::sweep_binder(beta, engine);
}
