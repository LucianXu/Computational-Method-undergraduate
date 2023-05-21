#include <bits/stdc++.h>

template <typename T>
struct mtx {
    int n, m;
    std::vector<std::vector<T>> a;

    mtx(int _n = 0, int _m = 0) {
        n = _n, m = _m;
        a = std::vector(n + 1, std::vector<T>(m + 1, (T) 0));
    }

    friend mtx matrixAdd(const mtx& A, const mtx& B, int n, int m) {
        mtx ans = A;
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                ans.a[i][j] += B.a[i][j];
            }
        }
        return ans;
    }

    friend mtx matrixSub(const mtx& A, const mtx& B, int n, int m) {
        mtx ans = A;
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                ans.a[i][j] -= B.a[i][j];
            }
        }
        return ans;
    }

    friend mtx matrixAdd(const mtx& A, const mtx& B, int n) { return matrixAdd(A, B, n, n); }

    friend mtx matrixSub(const mtx& A, const mtx& B, int n) { return matrixSub(A, B, n, n); }

    friend mtx matrixMul(const mtx& A, const T k, int n, int m) {
        mtx ans = A;
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                ans.a[i][j] *= k;
            }
        }
        return ans;
    }

    friend mtx matrixMul(const mtx& A, const T k, int n) { return matrixMul(A, k, n, n); }

    friend mtx matrixMul(const mtx& A, const mtx& B, int n, int m, int k) {
        mtx ans = mtx(n, k);
        for (int i = 1; i <= n; i++) {
            for (int l = 1; l <= m; l++) {
                T x = A.a[i][l];
                for (int j = 1; j <= k; j++) {
                    ans.a[i][j] += x * B.a[l][j];
                }
            }
        }
        return ans;
    }

    friend mtx matrixMul(const mtx& A, const mtx& B, int n) { return matrixMul(A, B, n, n, n); }

    friend std::vector<T> matrixVectorMul(const mtx& A, const std::vector<T>& x, int n, int m) {
        std::vector<T> y(m + 1);
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                y[i] += A.a[i][j] * x[j];
            }
        }
        return y;
    }

    friend std::vector<T> matrixVectorMul(const mtx& A, const std::vector<T>& x, int n) {
        return matrixVectorMul(A, x, n, n);
    }

    friend void subMatrix(mtx& A, const mtx& B, int n, int m, int le, int up) {
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                A.a[i][j] = B.a[up + i - 1][le + j - 1];
            }
        }
    }

    friend void lowerTriangularMatrixInv(mtx& A, int le, int up, int len) {
        /*
        | A   0 |       |        A^{-1}          0     |
        |       |  -->  |                              |
        | C   B |       |  -B^{-1} C A^{-1}    B^{-1}  |
        */
        if (len == 1) {
            A.a[le][up] = (T) 1 / A.a[le][up];
            return;
        }

        int l1 = len / 2, l2 = (len + 1) / 2;
        lowerTriangularMatrixInv(A, le, up, l1);
        lowerTriangularMatrixInv(A, le + l1, up + l1, l2);

        mtx AInv = mtx(l1, l1);
        mtx BInv = mtx(l2, l2);
        mtx C = mtx(l2, l1);

        subMatrix(AInv, A, l1, l1, le, up);
        subMatrix(BInv, A, l2, l2, le + l1, up + l1);
        subMatrix(C, A, l2, l1, le, up + l1);

        C = matrixMul(BInv, C, l2, l2, l1);
        C = matrixMul(C, AInv, l2, l1, l1);

        for (int i = 1; i <= l2; i++) {
            for (int j = 1; j <= l1; j++) {
                A.a[up + l1 + i - 1][le + j - 1] = -C.a[i][j];
            }
        }
    }

    friend mtx lowerTriangularMatrixInv(const mtx& A, int n) {
        mtx B = A;
        lowerTriangularMatrixInv(B, 1, 1, n);
        return B;
    }

    friend mtx diagonalMatrixInv(const mtx& A, int n) {
        mtx B = A;
        for (int i = 1; i <= n; i++) B.a[i][i] = (T) 1. / A.a[i][i];
        return B;
    }
};

typedef mtx<double> matrix;

std::vector<double> vectorAdd(const std::vector<double> A, const std::vector<double> B, int n) {
    std::vector<double> ans = A;
    for (int i = 1; i <= n; i++) ans[i] += B[i];
    return ans;
}

void Jacoib(
    const matrix& A, const std::vector<double>& b, const std::vector<double>& initVal, int n,
    int t) {
    matrix L = matrix(n, n);
    matrix D = matrix(n, n);
    matrix U = matrix(n, n);
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i > j) {
                L.a[i][j] = A.a[i][j];
            } else if (i == j) {
                D.a[i][j] = A.a[i][j];
            } else {
                U.a[i][j] = A.a[i][j];
            }
        }
    }

    matrix DInv = diagonalMatrixInv(D, n);
    matrix J = matrixMul(matrixMul(DInv, matrixAdd(L, U, n), n), -1., n);
    std::vector<double> f = matrixVectorMul(DInv, b, n);

    auto iterate = [&](const matrix& J, const std::vector<double>& f,
                       const std::vector<double>& x) -> std::vector<double> {
        return vectorAdd(matrixVectorMul(J, x, n), f, n);
    };

    std::vector<double> x = initVal;

    for (int _ = 1; _ <= t; _++) {
        x = iterate(J, f, x);
        std::cout << "This is the result of the " << _ << "-th iteration." << '\n';
        for (int i = 1; i <= n; i++) {
            std::cout << std::fixed << std::setprecision(8) << "x_" << i << " = " << x[i] << '\n';
        }
        std::cout << '\n';
    }
}

void Gauss_Seidel(
    const matrix& A, const std::vector<double>& b, const std::vector<double>& initVal, int n,
    int t) {
    matrix L = matrix(n, n);
    matrix D = matrix(n, n);
    matrix U = matrix(n, n);
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i > j) {
                L.a[i][j] = A.a[i][j];
            } else if (i == j) {
                D.a[i][j] = A.a[i][j];
            } else {
                U.a[i][j] = A.a[i][j];
            }
        }
    }

    matrix Inv = lowerTriangularMatrixInv(matrixAdd(D, L, n), n);
    matrix G = matrixMul(matrixMul(Inv, U, n), -1., n);
    std::vector<double> f = matrixVectorMul(Inv, b, n);

    auto iterate = [&](const matrix& G, const std::vector<double>& f,
                       const std::vector<double>& x) -> std::vector<double> {
        return vectorAdd(matrixVectorMul(G, x, n), f, n);
    };

    std::vector<double> x = initVal;

    std::cout << "This is the initial value." << '\n';
    for (int i = 1; i <= n; i++) {
        std::cout << std::fixed << std::setprecision(8) << "x_" << i << " = " << x[i] << '\n';
    }
    std::cout << '\n';

    for (int _ = 1; _ <= t; _++) {
        x = iterate(G, f, x);
        std::cout << "This is the result of the " << _ << "-th iteration." << '\n';
        for (int i = 1; i <= n; i++) {
            std::cout << std::fixed << std::setprecision(8) << "x_" << i << " = " << x[i] << '\n';
        }
        std::cout << '\n';
    }
}

int main() {
    // int n;
    // std::cin >> n;
    // matrix A = matrix(n, n);
    // for (int i = 1; i <= n; i++) {
    //     for (int j = 1; j <= n; j++) {
    //         std::cin >> A.a[i][j];
    //     }
    // }
    // std::vector<double> b(n + 1), x(n + 1);
    // for (int i = 1; i <= n; i++) std::cin >> b[i];

    int n;
    std::cin >> n;
    matrix A = matrix(n, n);
    for (int i = 1; i <= n; i++) {
        A.a[i][i] = 0.5;
        if (i + 1 <= n) A.a[i][i + 1] = A.a[i + 1][i] = 0.5;
        if (i + 2 <= n) A.a[i][i + 2] = A.a[i + 2][i] = 0.5;
        if (i <= n / 2) A.a[i][2 * i] = A.a[2 * i][i] = 0.5;
    }
    std::vector<double> b(n + 1), x(n + 1);
    for (int i = 1; i <= n; i++) b[i] = 1;

    // Jacoib(A, b, x, n, 5);
    Gauss_Seidel(A, b, x, n, 5);

    return 0;
}
