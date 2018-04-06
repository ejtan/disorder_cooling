#include <vector>
#include <algorithm>

#include "utility.h"


template <class Model>
std::vector<double> compute_binder(const std::vector<double> &T, Model &model)
{
    std::vector<double> binder(T.size());
    std::fill(binder.begin(), binder.end(), 0);
    run_mc_binder(T, binder, model);

    return binder;
}


template <class Model>
std::vector<double> compute_binder(const std::vector<double> &T, Model &model,
        double delta, int n_run)
{
    std::vector<double> binder(T.size());
    std::fill(binder.begin(), binder.end(), 0);

    for (int run = 0; run < n_run; run++) {
        model.set_exchange(delta);
        run_mc_binder(T, binder, model);
    } // Loop over runs.

    // Normalize data
    std::transform(binder.begin(), binder.end(), binder.begin(),
            [n_run](double val)->double { return val / static_cast<double>(n_run); });

    return binder;
}
