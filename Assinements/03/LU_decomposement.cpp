#include <bits/stdc++.h>

template <typename T>
struct mtx {
    int n, m;
    std::vector<std::vector<T>> a;

    mtx(int _n = 0, int _m = 0) {
        n = _n, m = _m;
        a = std::vector(n + 1, std::vector<T>(m + 1, (T) 0));
    }

    friend void subMatrix(mtx& ma, const mtx& mb, int n, int m, int le, int up) {
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                ma.a[i][j] = mb.a[up + i - 1][le + j - 1];
            }
        }
    }

    friend mtx matrixMul(const mtx& a, const mtx& b, int n) {
        mtx c = mtx(n, n);
        for (int i = 1; i <= n; i++) {
            for (int l = 1; l <= n; l++) {
                T x = a[i][l];
                for (int j = 1; j <= n; j++) {
                    c[i][j] += x * b[l][j];
                }
            }
        }
        return c;
    }

    friend mtx matrixMul(const mtx& ma, const mtx& mb, int n, int m, int k) {
        mtx c = mtx(n, k);
        for (int i = 1; i <= n; i++) {
            for (int l = 1; l <= m; l++) {
                T x = ma.a[i][l];
                for (int j = 1; j <= k; j++) {
                    c.a[i][j] += x * mb.a[l][j];
                }
            }
        }
        return c;
    }

    friend void lowerTriangularMatrixInv(mtx& ma, int le, int up, int len) {
        /*
        | A   0 |       |        A^{-1}          0     |
        |       |  -->  |                              |
        | C   B |       |  -B^{-1} C A^{-1}    B^{-1}  |
        */
        if (len == 1) {
            ma.a[le][up] = (T) 1 / ma.a[le][up];
            return;
        }

        int l1 = len / 2, l2 = (len + 1) / 2;
        lowerTriangularMatrixInv(ma, le, up, l1);
        lowerTriangularMatrixInv(ma, le + l1, up + l1, l2);

        mtx AInv = mtx(l1, l1);
        mtx BInv = mtx(l2, l2);
        mtx C = mtx(l2, l1);

        subMatrix(AInv, ma, l1, l1, le, up);
        subMatrix(BInv, ma, l2, l2, le + l1, up + l1);
        subMatrix(C, ma, l2, l1, le, up + l1);

        C = matrixMul(BInv, C, l2, l2, l1);
        C = matrixMul(C, AInv, l2, l1, l1);

        for (int i = 1; i <= l2; i++) {
            for (int j = 1; j <= l1; j++) {
                ma.a[up + l1 + i - 1][le + j - 1] = -C.a[i][j];
            }
        }
    }
    friend std::pair<mtx, mtx> tridiagonalMatrix_LU_Decomposement(const mtx& A, int n) {
        mtx L = mtx(n, n);
        mtx U = A;
        for (int i = 1; i <= n; i++) L.a[i][i] = 1;
        for (int i = 1; i < n; i++) {
            T inv = U.a[i + 1][i] / U.a[i][i];
            for (int j = 0; j <= 1; j++) {
                U.a[i + 1][i + j] -= U.a[i][i + j] * inv;
            }
            for (int j = 1; j <= i; j++) {
                L.a[i + 1][j] -= L.a[i][j] * inv;
            }
        }
        // lowerTriangularMatrixInv(L, 1, 1, n);
        return std::make_pair(L, U);
    }

    friend std::vector<T> tridiagonalMatrixSolve(const mtx& A, const std::vector<T>& b, int n) {
        auto [LInv, U] = tridiagonalMatrix_LU_Decomposement(A, n);
        std::vector<T> bb(n + 1);
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                bb[i] += LInv.a[i][j] * b[j];
            }
        }
        std::vector<T> x(n + 1);
        for (int i = n; i >= 1; i--) {
            x[i] = bb[i];
            for (int j = i + 1; j <= n; j++) {
                x[i] -= x[j] * U.a[i][j];
            }
            x[i] /= U.a[i][i];
        }
        return x;
    }

    friend std::pair<mtx, mtx> cycleTridiagonalMatrix_LU_Decomposement(const mtx& A, int n) {
        mtx L = mtx(n, n);
        mtx U = A;
        for (int i = 1; i <= n; i++) L.a[i][i] = 1;
        for (int i = 1; i < n; i++) {
            if (i != n - 1) {
                T inv1 = U.a[i + 1][i] / U.a[i][i];
                T inv2 = U.a[n][i] / U.a[i][i];
                for (int j = 0; j <= 1; j++) {
                    U.a[i + 1][i + j] -= U.a[i][i + j] * inv1;
                    U.a[n][i + j] -= U.a[i][i + j] * inv2;
                }
                for (int j = 1; j <= i; j++) {
                    L.a[i + 1][j] -= L.a[i][j] * inv1;
                    L.a[n][j] -= L.a[i][j] * inv2;
                }
            } else {
                T inv = U.a[i + 1][i] / U.a[i][i];
                for (int j = 0; j <= 1; j++) {
                    U.a[i + 1][i + j] -= U.a[i][i + j] * inv;
                }
                for (int j = 1; j <= i; j++) {
                    L.a[i + 1][j] -= L.a[i][j] * inv;
                }
            }
        }
        // lowerTriangularMatrixInv(L, 1, 1, n);
        return std::make_pair(L, U);
    }

    friend std::vector<T> cycleTridiagonalMatrixSolve(
        const mtx& A, const std::vector<T>& b, int n) {
        auto [LInv, U] = cycleTridiagonalMatrix_LU_Decomposement(A, n);
        std::vector<T> bb(n + 1);
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                bb[i] += LInv.a[i][j] * b[j];
            }
        }
        std::vector<T> x(n + 1);
        for (int i = n; i >= 1; i--) {
            x[i] = bb[i];
            for (int j = i + 1; j <= n; j++) {
                x[i] -= x[j] * U.a[i][j];
            }
            x[i] /= U.a[i][i];
        }
        return x;
    }
};

typedef mtx<double> matrix;

int main() {
    // First, we test tridiagona matrix. //

    int n;
    std::cin >> n;
    matrix A = matrix(n, n);
    for (int i = 1; i <= n; i++) std::cin >> A.a[i][i];
    for (int i = 1; i <= n - 1; i++) std::cin >> A.a[i][i + 1];
    for (int i = 2; i <= n; i++) std::cin >> A.a[i][i - 1];
    std::vector<double> b(n + 1);
    for (int i = 1; i <= n; i++) std::cin >> b[i];

    auto X = tridiagonalMatrixSolve(A, b, n);
    for (int i = 1; i <= n; i++) {
        std::cout << std::fixed << std::setprecision(6) << "x_i = " << X[i] << '\n';
    }


    // Then, we test cycle tridiagona matrix. //
    for (int i = 1; i <= n; i++) std::cin >> A.a[i][i];
    for (int i = 1; i <= n; i++) std::cin >> A.a[i][i < n ? i + 1 : 1];
    for (int i = 1; i <= n; i++) std::cin >> A.a[i < n ? i + 1 : 1][i];
    for (int i = 1; i <= n; i++) std::cin >> b[i];

    X = cycleTridiagonalMatrixSolve(A, b, n);
    for (int i = 1; i <= n; i++) {
        std::cout << std::fixed << std::setprecision(6) << "x_i = " << X[i] << '\n';
    }

    return 0;
}
