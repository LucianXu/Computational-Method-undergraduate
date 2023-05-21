import numpy as np


class Matrix:
    def __init__(self, n=0, m=0):
        self.n = n
        self.m = m
        self.a = np.zeros([n, m])

    def tridiagonalMatrix_LU_Decomposement(self, n):
        L = Matrix(n, n)
        U = Matrix(n, n)
        U.a = self.a[:]
        for i in range(0, n):
            L.a[i][i] = 1
        for i in range(0, n - 1):
            inv = U.a[i + 1][i] / U.a[i][i]
            for j in range(0, 2):
                U.a[i + 1][i + j] -= U.a[i][i + j] * inv
            for j in range(0, i + 1):
                L.a[i + 1][j] -= L.a[i][j] * inv
        return L, U

    def tridiagonalMatrixSolve(self, b, n):
        LInv, U = self.tridiagonalMatrix_LU_Decomposement(n)
        b = np.matmul(LInv.a, b)
        x = np.zeros([n, 1])
        for i in range(n - 1, -1, -1):
            x[i, 0] = b[i, 0]
            for j in range(i + 1, n):
                x[i, 0] -= x[j, 0] * U.a[i, j]
            x[i, 0] /= U.a[i, i]
        return x

    def cycleTridiagonalMatrix_LU_Decomposement(self, n):
        L = Matrix(n, n)
        U = Matrix(n, n)
        U.a = self.a[:]
        for i in range(0, n):
            L.a[i][i] = 1
        for i in range(0, n - 1):
            if i != n - 2:
                inv1 = U.a[i + 1][i] / U.a[i][i]
                inv2 = U.a[n - 1][i] / U.a[i][i]
                for j in range(0, 2):
                    U.a[i + 1][i + j] -= U.a[i][i + j] * inv1
                    U.a[n - 1][i + j] -= U.a[i][i + j] * inv2
                for j in range(0, i + 1):
                    L.a[i + 1][j] -= L.a[i][j] * inv1
                    L.a[n - 1][j] -= L.a[i][j] * inv2
            else:
                inv = U.a[i + 1][i] / U.a[i][i]
                for j in range(0, 2):
                    U.a[i + 1][i + j] -= U.a[i][i + j] * inv
                for j in range(0, i + 1):
                    L.a[i + 1][j] -= L.a[i][j] * inv
        return L, U

    def cycleTridiagonalMatrixSolve(self, b, n):
        LInv, U = self.cycleTridiagonalMatrix_LU_Decomposement(n)
        b = np.matmul(LInv.a, b)
        x = np.zeros([n, 1])
        for i in range(n - 1, -1, -1):
            x[i, 0] = b[i, 0]
            for j in range(i + 1, n):
                x[i, 0] -= x[j, 0] * U.a[i, j]
            x[i, 0] /= U.a[i, i]
        return x


# First, we test tridiagona matrix.
n = int(input())
A = Matrix(n, n)

ipt = input().split(" ")
for i in range(0, n):
    A.a[i][i] = ipt[i]

ipt = input().split(" ")
for i in range(0, n - 1):
    A.a[i][i + 1] = ipt[i]

ipt = input().split(" ")
for i in range(0, n - 1):
    A.a[i + 1][i] = ipt[i]

b = np.zeros([n, 1])
ipt = input().split(" ")
for i in range(0, n):
    b[i, 0] = ipt[i]

x = A.tridiagonalMatrixSolve(b, n)
for i in range(0, n):
    print("x_{} = {}".format(i + 1, x[i, 0]))

# Then, we test cycle tridiagona matrix.
n = int(input())
A = Matrix(n, n)

ipt = input().split(" ")
for i in range(0, n):
    A.a[i][i] = ipt[i]

ipt = input().split(" ")
for i in range(0, n):
    A.a[i][(i + 1) % n] = ipt[i]

ipt = input().split(" ")
for i in range(0, n):
    A.a[(i + 1) % n][i] = ipt[i]

b = np.zeros([n, 1])
ipt = input().split(" ")
for i in range(0, n):
    b[i, 0] = ipt[i]

x = A.cycleTridiagonalMatrixSolve(b, n)
for i in range(0, n):
    print("x_{} = {}".format(i + 1, x[i, 0]))
