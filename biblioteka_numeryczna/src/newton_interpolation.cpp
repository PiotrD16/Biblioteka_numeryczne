#include "../include/newton_interpolation.h"

Interpolation::Interpolation(const std::vector<double>& x_vals, const std::vector<double>& fx_vals)
    : x_values(x_vals), fx_values(fx_vals) {
    computeNewtonCoefficients();
}

void Interpolation::computeNewtonCoefficients() {
    int n = x_values.size();
    std::vector<std::vector<double>> dd(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        dd[i][0] = fx_values[i];
    }

    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            dd[i][j] = (dd[i + 1][j - 1] - dd[i][j - 1]) / (x_values[i + j] - x_values[i]);
        }
    }

    coefficients.resize(n);
    for (int i = 0; i < n; i++) {
        coefficients[i] = dd[0][i];
    }
}

double Interpolation::evaluateNewton(double x) const {
    double result = coefficients.back();
    for (int i = coefficients.size() - 2; i >= 0; --i) {
        result = result * (x - x_values[i]) + coefficients[i];
    }
    return result;
}

std::vector<double> Interpolation::getCoefficients() const {
    return coefficients;
}