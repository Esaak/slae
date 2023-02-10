//
// Created by ilya on 05.02.23.
//
#include "tridiagonal_matrix.hpp"
#ifndef SLAE_SLAE_SOLVER_HPP
#define SLAE_SLAE_SOLVER_HPP

template <typename T>
std::vector<T> tridiagonal_matrix_algorithm(Tridiagonal_matrix<T> A, std::vector<T> D){
    int N = static_cast<int>(D.size()) - 1;
    A(0).c/=A(0).b;
    D[0]/=A(0).b;
    for(int i = 1; i < N; i++){
        A(i).c /= A(i).b - A(i).a * A(i-1).c;
        D[i] = (D[i] - A(i).a * D[i-1]) / (A(i).b - A(i).a * A(i-1).c);
    }
    D.back() = (D.back() - A(N).a * D[N - 1]) / (A(N).b - A(N).a * A(N-1).c);

    for(int i = N; i-- > 0; ){
        D[i] -=A(i).c * D[i+1];
    }
    return D;
}

#endif //SLAE_SLAE_SOLVER_HPP
