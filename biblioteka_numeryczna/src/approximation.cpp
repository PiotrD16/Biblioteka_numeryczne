#include "../include/approximation.h"

double approximation::monomial(int degree, double x) 
{
    return std::pow(x, degree);
}

double approximation::integrate(std::function<double(double)> func, int n = 1000) {

    if (n % 2 == 1) ++n; // musi być parzyste
    double h = (b - a) / n;
    double sum = func(a) + func(b);

    for (int i = 1; i < n; i += 2) {
        sum += 4 * func(a + i * h);
    }
    for (int i = 2; i < n; i += 2) {
        sum += 2 * func(a + i * h);
    }

    return (h / 3.0) * sum;
}

void approximation::least_squares_approximation(std::function<double(double)> f) {
    std::vector<std::vector<double>> A(degree + 1, std::vector<double>(degree + 1));
    std::vector<double> b_vec(degree + 1);

    // Wypełnianie macierzy A i wektora b
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            A[i][j] = integrate([=](double x) {
                return monomial(i, x) * monomial(j, x);
            });
        }
        b_vec[i] = integrate([=](double x) {
            return f(x) * monomial(i, x);
        });
    }

    // Eliminacja Gaussa – eliminacja do postaci górnotrójkątnej
    for (int i = 0; i <= degree; ++i) {
        double diag = A[i][i];
        for (int j = 0; j <= degree; ++j) A[i][j] /= diag;
        b_vec[i] /= diag;

        for (int k = i + 1; k <= degree; ++k) {
            double factor = A[k][i];
            for (int j = 0; j <= degree; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b_vec[k] -= factor * b_vec[i];
        }
    }

    // Podstawianie wstecz
    std::vector<double> c(degree + 1);
    for (int i = degree; i >= 0; --i) {
        c[i] = b_vec[i];
        for (int j = i + 1; j <= degree; ++j) {
            c[i] -= A[i][j] * c[j];
        }
    }

    // Wypisanie współczynników
    std::cout << "\nWspółczynniki aproksymacji:\n";
    for (int i = 0; i <= degree; ++i) {
        std::cout << "c[" << i << "] = " << c[i] << "\n";
    }
}
