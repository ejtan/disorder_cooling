#include <vector>
#include <random>
#include <omp.h>


template <class Model>
void run_mc_energy(const std::vector<double> &T, std::vector<double> &E, Model model)
{
    #pragma omp parallel firstprivate(model) num_threads(4)
    {
        int chunk = T.size() / 4;
        int thd_id = omp_get_thread_num();
        std::random_device rd;
        std::mt19937 engine(rd());

        for (int i = chunk * thd_id; i < (chunk * (thd_id + 1)); i++) {
            model.set_spin();
            E[i] += model.sweep_energy(1.0 / T[i], engine);
        } // Loop over chunk
    } // Parallel region
}


template <class Model>
void run_mc_binder(const std::vector<double> &T, std::vector<double> &binder, Model model)
{
    #pragma omp parallel firstprivate(model) num_threads(4)
    {
        int chunk = T.size() / 4;
        int thd_id = omp_get_thread_num();
        std::random_device rd;
        std::mt19937 engine(rd());

        for (int i = chunk * thd_id; i < (chunk * (thd_id + 1)); i++) {
            model.set_spin();
            binder[i] += model.sweep_binder(1.0 / T[i], engine);
        } // Loop over chunk
    } // Parallel region
}


double trapezoid(const std::vector<double> &x, const std::vector<double> &y, int idx)
{
    int N = x.size();
    double sum = 0.0;

    #pragma omp simd reduction(+:sum)
    for (int i = idx; i < N - 1; i++) {
        double dx = x[i + 1] - x[i];
        sum += dx * (y[i] + y[i + 1]);
    } // Loop to compute integral

    return 0.5 * sum;
}

