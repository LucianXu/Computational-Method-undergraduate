#include <bits/stdc++.h>

#include "Matrix.h"

typedef mtx<double> matrix;


std::vector<std::array<double, 4>> cubicSplineInterpolation(
    const std::vector<double>& X, const std::vector<double>& Y, int n, double bound_1,
    double bound_2) {
    std::vector<std::array<double, 4>> answer;

    std::vector<double> mu(n + 1), lambda(n + 1);
    for (int i = 2; i <= n; i++) {
        mu[i] = X[i] - X[i - 1] / (X[i + 1] - X[i - 1]);
        lambda[i] = X[i + 1] - X[i] / (X[i + 1] - X[i - 1]);
    }

    auto secondOrderDifferece = [&](const int& x, const int& y, const int& z) -> double {
        double f1 = (Y[y] - Y[x]) / (X[y] - X[x]);
        double f2 = (Y[z] - Y[y]) / (X[z] - X[y]);
        return (f2 - f1) / (X[z] - X[x]);
    };

    std::vector<double> d(n + 2);
    for (int i = 1; i <= n; i++) {
        if (i == 1) {
            double f = (Y[2] - Y[1]) / (X[2] - X[1]);
            d[i] = (f - bound_1) / (X[2] - X[1]);
        } else if (i == n + 1) {
            double f = (Y[n + 1] - Y[n]) / (X[n + 1] - X[n]);
            d[i] = (bound_2 - f) / (X[n + 1] - X[n]);
        } else {
            d[i] = secondOrderDifferece(i - 1, i, i + 1);
        }
    }

    matrix A = matrix(n + 1, n + 1);
    for (int i = 1; i <= n + 1; i++) {
        A.a[i][i] = 2.;
        if (i == 1) {
            A.a[i][i + 1] = 1.;
        } else if (i == n + 1) {
            A.a[i][i - 1] = 1.;
        } else {
            A.a[i][i + 1] = lambda[i - 1];
            A.a[i][i - 1] = mu[i - 1];
        }
    }

    auto M = tridiagonalMatrixSolve(A, d, n + 1);

    for (int i = 1; i <= n; i++) {
        std::array<double, 4> func;
        func[0] = Y[i];
        func[1] = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]) - (X[i + 1] - X[i]) / 2 * M[i] -
                  (X[i + 1] - X[i]) / 6 * (M[i + 1] - M[i]);
        func[2] = M[i] / 2;
        func[3] = (M[i + 1] - M[i]) / (6 * (X[i + 1] - X[i]));
        answer.push_back(func);
    }
    return answer;
}

int main() {
    // input n + 1 pairs of value, from 1 to n + 1 //

    int n;
    std::cin >> n;
    std::vector<double> X(n + 2), Y(n + 2);
    double bound_1, bound_2;
    for (int i = 1; i <= n + 1; i++) std::cin >> X[i];
    for (int i = 1; i <= n + 1; i++) std::cin >> Y[i];
    std::cin >> bound_1 >> bound_2;

    auto answer = cubicSplineInterpolation(X, Y, n, bound_1, bound_2);
    for (const auto& func : answer) {
        std::cerr << func[0] << ' ' << func[1] << ' ' << func[2] << ' ' << func[3] << '\n';
    }
    return 0;
}
