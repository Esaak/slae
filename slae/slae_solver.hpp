//
// Created by ilya on 05.02.23.
//
#include "tridiagonal_matrix.hpp"
#ifndef SLAE_SLAE_SOLVER_HPP
#define SLAE_SLAE_SOLVER_HPP

template <typename T>
std::vector<T> tridiagonal_matrix_algorithm(Tridiagonal_matrix<T>& A, std::vector<T> D){
    int N = static_cast<int>(D.size()) - 1;
    std::vector<T> a(N+1);
    std::vector<T> b(N+1);
    std::vector<T> c(N+1);
    for(int i = 0; i <= N; i++){
        a[i] = A(i).a;
        b[i] = A(i).b;
        c[i] = A(i).c;
    }
    c[0]/=b[0];
    D[0]/=b[0];
    for(int i = 1; i < N; i++){
        c[i] /= b[i] - a[i] * c[i-1];
        D[i] = (D[i] - a[i] * D[i-1]) / (b[i] - a[i] * c[i-1]);
    }
    D.back() = (D.back() - a[N] * D[N - 1]) / (b[N] - a[N] * c[N-1]);

    for(int i = N; i-- > 0; ){
        D[i] -=c[i] * D[i+1];
    }
    return D;
}


#endif //SLAE_SLAE_SOLVER_HPP
