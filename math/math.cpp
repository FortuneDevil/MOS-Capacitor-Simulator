#include "math.hpp"

// Bernoulli function = x/(e^x-1)
double B(double x){
    return bernoulli(x);
}

// Derivative of the Bernoulli function
double dB(double x){
    return dBernoulli(x);
}

std::vector<double> solveFluxSystem(const std::vector<double>& up, const std::vector<double>& down,
                                    const std::vector<double>& sink, const std::vector<double>& source, double last_value){
    int n = sink.size();
    std::vector<double> y(n), e(n), beta(n);

    // Forward Elimination. After it, node k obeys: up[k] * y[k+1] - (down[k] + e[k]) * y[k] = -beta[k]
    e[0] = sink[0];
    beta[0] = source[0];
    for (int k = 1; k < n - 1; k++){
        double pivot = down[k - 1] + e[k - 1];
        e[k] = sink[k] + up[k - 1] * e[k - 1] / pivot;
        beta[k] = source[k] + down[k - 1] * beta[k - 1] / pivot;
    }

    // Backward Substitution
    y[n - 1] = last_value;
    for (int k = n - 2; k >= 0; k--){
        y[k] = (up[k] * y[k + 1] + beta[k]) / (down[k] + e[k]);
    }
    return y;
}
