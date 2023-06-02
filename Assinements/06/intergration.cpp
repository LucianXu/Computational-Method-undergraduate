#include <bits/stdc++.h>

double f(double x) {
    // todo //
    return std::sin(x) / x;
}

double intergration(std::pair<double, double> interval, int n, int type) {
    auto [a, b] = interval;
    double h = (b - a) / n;
    double res = 0;
    if (type <= 3) {
        for (int i = 0; i < n; i++) {
            double le = a + i * h;
            double ri = a + (i + 1) * h;
            if (type == 1) {    // trapezoid //
                res += h * (f(le) + f(ri)) / 2;
            } else if (type == 2) {    // Simpson //
                res += h * (f(le) + 4 * f((le + ri) / 2) + f(ri)) / 6;
            } else if (type == 3) {    // Cotes //
                res += h *
                       (7 * f(le) + 32 * f((3 * le + ri) / 4) + 12 * f((le + ri) / 2) +
                        32 * f((le + 3 * ri) / 4) + 7 * f(ri)) /
                       90;
            }
        }
    } else if (type == 4) {    // Gauss-Legender //
        auto getLegenderValue = [&](double x, int n) -> double {
            if (n == 0) {
                return 1.0;
            } else if (n == 1) {
                return x;
            }
            std::vector<double> value(n + 1);
            value[0] = 1, value[1] = x;
            for (int i = 2; i <= n; i++) {
                value[i] = x * (2 * i - 1) * value[i - 1] - (i - 1) * value[i - 2];
                value[i] /= i;
            }
            return value.back();
        };
        auto getZeros = [&](int n) -> std::vector<double> {
            std::vector<double> roots = {-1, 0, 1}, roots_tmp;
            for (int i = 2; i <= n; i++) {
                roots_tmp = {-1};
                for (int j = 0; j < (int) roots.size() - 1; j++) {
                    double l = roots[j], r = roots[j + 1];
                    while ((r - l) > 1e-10) {
                        double mid = (l + r) / 2;
                        double l_val = getLegenderValue(l, i);
                        double mid_val = getLegenderValue(mid, i);
                        if (fabs(mid) <= 1e-10) {
                            l = mid;
                            break;
                        } else if ((l_val < 0 and mid_val < 0) or (l_val > 0 and mid_val > 0)) {
                            l = mid;
                        } else {
                            r = mid;
                        }
                    }
                    roots_tmp.push_back(l);
                }
                roots_tmp.push_back(1.0);
                for (auto& x : roots_tmp) {
                    if (fabs(x) < 1e-10) x = 0;
                }
                std::swap(roots, roots_tmp);
            }
            roots.pop_back();
            reverse(roots.begin(), roots.end());
            roots.pop_back();
            reverse(roots.begin(), roots.end());
            return roots;
        };
        auto getLegenderDerivationValue = [&](double x, int n) -> double {
            if (n == 0) {
                return 0.0;
            } else if (n == 1) {
                return 1.0;
            }
            std::vector<double> value(n + 1), derivation(n + 1);
            value[0] = 1, value[1] = x;
            derivation[0] = 0.0, derivation[1] = 1.0;
            for (int i = 2; i <= n; i++) {
                value[i] = x * (2 * i - 1) * value[i - 1] - (i - 1) * value[i - 2];
                value[i] /= i;
                derivation[i] = x * (2 * i - 1) * derivation[i - 1] - (i - 1) * derivation[i - 2] +
                                (2 * i - 1) * value[i - 1];
                derivation[i] /= i;
            }
            return derivation.back();
        };
        auto roots = getZeros(n + 1);
        std::vector<double> weight;
        for (const auto& x : roots) {
            double val = getLegenderDerivationValue(x, n + 1);
            val = 2 / (1 - x * x) / (val * val);
            weight.push_back(val);
        }
        for (int i = 0; i <= n; i++) {
            double x = (a + b) / 2 + (b - a) / 2 * roots[i];
            res += (b - a) / 2 * weight[i] * f(x);
        }
    }
    return res;
};


int main() {
    std::cin.tie(nullptr)->sync_with_stdio(false);
    std::cout << std::fixed << std::setprecision(14);

    std::cout << "For the function sin(x) / x and the integral interval [1, 5], when the "
                 "interval is divided into n sections, "
                 "here n takes values from 3 to 7, the four methods of calculating the integral "
                 "give the following results: "
              << '\n';
    for (int i = 3; i <= 7; i++) {
        std::cout << "When [3, 5] is divided into " << i << " sections, " << '\n';
        std::cout << "here is the result of trapezoid formular: " << intergration({1, 5}, i, 1)
                  << ",\n";
        std::cout << "here is the result of Simpson formular: " << intergration({1, 5}, i, 2)
                  << ",\n";
        std::cout << "here is the result of Cotes formular: " << intergration({1, 5}, i, 3)
                  << ",\n";
        std::cout << "here is the result of Gauss-Legender formular: " << intergration({1, 5}, i, 4)
                  << ",\n";

        std::cout << '\n';
    }
    return 0;
}
