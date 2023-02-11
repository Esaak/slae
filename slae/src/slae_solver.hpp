//
// Created by ilya on 05.02.23.
//
#ifndef SLAE_SLAE_SOLVER_HPP
#define SLAE_SLAE_SOLVER_HPP
#include "tridiagonal_matrix.hpp"

namespace solvers {
    template<typename T>
    std::vector<T> tridiagonal_matrix_algorithm(const Tridiagonal_matrix<T> &A, const std::vector<T> &D) {
        std::vector<T> X = D;
        std::vector<T> Temp(D.size());

        int N = static_cast<int>(D.size()) - 1;

        Temp[0] = A(0, 0 + 1) / A(0, 0);
        X[0] /= A(0, 0);
        for (int i = 1; i < N; i++) {
            Temp[i] = A(i, i + 1) /( A(i, i) - A(i, i - 1) * Temp[i - 1]);
            X[i] = (X[i] - A(i, i - 1) * X[i - 1]) / (A(i, i) - A(i, i - 1) * Temp[i - 1]);
        }
        X.back() = (X.back() - A(N, N - 1) * X[N - 1]) / (A(N, N) - A(N, N - 1) * Temp[N - 1]);

        for (int i = N; i-- > 0;) {
            X[i] -= Temp[i] * X[i + 1];
        }
        return X;
    }
}
#endif //SLAE_SLAE_SOLVER_HPP
