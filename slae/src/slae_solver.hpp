//
// Created by ilya on 05.02.23.
//
#ifndef SLAE_SLAE_SOLVER_HPP
#define SLAE_SLAE_SOLVER_HPP

#include "tridiagonal_matrix.hpp"
#include "CSR_matrix.hpp"
namespace solvers {
    template<std::floating_point T>
    std::vector<T> tridiagonal_matrix_algorithm(const Tridiagonal_matrix<T> &A, const std::vector<T> &D) {
        std::vector<T> X = D;
        std::vector<T> Temp(D.size());

        int N = static_cast<int>(D.size()) - 1;

        Temp[0] = A(0, 0 + 1) / A(0, 0);
        X[0] /= A(0, 0);
        for (int i = 1; i < N; i++) {
            Temp[i] = A(i, i + 1) / (A(i, i) - A(i, i - 1) * Temp[i - 1]);
            X[i] = (X[i] - A(i, i - 1) * X[i - 1]) / (A(i, i) - A(i, i - 1) * Temp[i - 1]);
        }
        X.back() = (X.back() - A(N, N - 1) * X[N - 1]) / (A(N, N) - A(N, N - 1) * Temp[N - 1]);

        for (int i = N; i-- > 0;) {
            X[i] -= Temp[i] * X[i + 1];
        }
        return X;
    }


    namespace _diny {
        template<typename T, IsArithmetical<T> = true>
        T sgn(T val) {
            return (T(0) < val) - (val < T(0));
        }

        template<std::floating_point T>
        std::vector<T> orthogonal_vector(const std::vector<T> &x) {
            T e;
            std::vector<T> new_x = x;
            e = sgn(new_x[0]) * sqrt(std::inner_product(new_x.begin(), new_x.end(), new_x.begin(), T(0)));
            new_x[0] += e;
            return new_x;
        }

        template<std::floating_point T>
        T search_projection(const std::vector<T> &v, const std::vector<T> &a, const T norm) {
            T betta = std::inner_product(v.begin(), v.end(), a.begin(), T(0));
            return betta / norm;
        }

        //this is experimental function does not use
        template<std::floating_point T>
        void find_Q_test(Mrx::Matrix<T> &Q, const std::vector<T> &v) {
            Mrx::Matrix<T> temp_Q = Q;
            for (std::size_t i = Q.get_column_size() - v.size(); i < Q.get_column_size(); i++) {
                T projection = search_projection(v, temp_Q.get_row(i, Q.get_column_size() - v.size(),
                                                                   Q.get_column_size()));
                for (std::size_t j = 0; j < Q.get_column_size(); j++) {
                    Q(i, j) -= 2 * v[j - (Q.get_column_size() - v.size())] * projection;
                }
            }
        }
    }

    //do not use please, for use this function you need save all ortoghonal vectors
    template<std::floating_point T>
    Mrx::Matrix<T> explicit_multiplication_for_Q(const std::vector<std::vector<T>> &matrix) {
        using namespace _diny;
        std::size_t N = matrix[0].size();
        Mrx::Matrix<T> temp;
        temp = temp.eye(N, N);
        for (std::size_t p = 0; p < matrix.size() - 1; p++) {
            Mrx::Matrix<T> temp1;
            Mrx::Matrix<T> temp2;
            temp1 = temp1.eye(N);
            temp2 = temp2.eye(matrix[p].size());
            std::size_t i;
            T norm = std::inner_product(matrix[p].begin(), matrix[p].end(), matrix[p].begin(), T(0));
            for (i = N - matrix[p].size(); i < N; ++i) {
                T projection = search_projection(matrix[p], temp2.get_column(i - (N - matrix[p].size()), 0,
                                                                             temp2.get_column_size()), norm);
                for (std::size_t j = N - matrix[p].size(); j < N; j++) {
                    temp1(j, i) -=
                            2 * matrix[p][j - (N - matrix[p].size())] * projection;//??????????????, ?????????? ???????????????? ???? ?? ????????
                }
            }
            temp *= temp1;
        }
        return temp;
    }

    template<std::floating_point T>
    std::pair<Mrx::Matrix<T>, Mrx::Matrix<T>> Householder(const Mrx::Matrix<T> &A) {
        using namespace _diny;
        Mrx::Matrix<T> new_matrix = A;
        Mrx::Matrix<T> Q;
        Q = Q.eye(A.get_column_size());
        for (std::size_t i = 0; i + 1 < new_matrix.get_column_size(); ++i) {
            std::vector<T> temp1 = new_matrix.get_column(i, i, new_matrix.get_column_size());
            std::vector<T> ort_vector = orthogonal_vector(temp1);
            T norma = std::inner_product(ort_vector.begin(), ort_vector.end(), ort_vector.begin(), T(0));
            for (std::size_t j = i; j < new_matrix.get_column_size(); j++) {
                if (i == j) {
                    T projection = T(0.5);
                    std::size_t p = j;
                    for (; p <= i; p++) {
                        new_matrix(p, j) -= 2 * ort_vector[p - j] * projection;
                    }
                    for (; p < new_matrix.get_column_size(); p++) {
                        new_matrix(p, j) = 0;
                    }
                } else {
                    T projection = search_projection(ort_vector, new_matrix.get_column(j, i, new_matrix.get_column_size()), norma);
                    for (std::size_t p = 0; p + i < new_matrix.get_column_size(); p++) {
                        new_matrix(p + i, j) -= 2 * ort_vector[p] * projection;
                    }
                }
            }
            for (std::size_t p = 0; p < Q.get_column_size(); p++) {
                T projection = search_projection(ort_vector, Q.get_column(p, i, Q.get_column_size()), norma);
                for (std::size_t j = 0; j + i < Q.get_column_size(); j++) {
                    Q(j + i, p) -= 2 * ort_vector[j] * projection;
                }
            }
        }
        Mrx::Matrix<T> QT = Q.transponse();
        return std::make_pair(QT, new_matrix);
    }
    /*
    template<std::floating_point T>
    T vector_norm(const std::vector<T>& vec){
        return std::inner_product(vec.begin(), vec.end(), vec.begin(), T(0));
    }

    template<std::floating_point T>
    std::vector<T>discrepancy(const CSR_matrix_space::CSR_matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x){
        //b- Ax;
        std::vector<T> ax = A * x;
        std::vector<T> answ(b.size());
        for(std::size_t i = 0; i < b.size(); i++){
            answ = b[i] - ax[i];
        }
        return answ;
    }*/
    /*
    template<std::floating_point T>
    std::vector<T>MPI (const CSR_matrix_space::CSR_matrix<T>& A, const std::vector<T>& b, T tau, T discrepancy0, const std::vector<T>& x0){
        std::vector<T> x = x0;
        std::vector<T> x_next = x0;
        while(x_next == x0 || discrepancy(A, b, x) >= discrepancy0 ){
            for(std::size_t i = 0; i < x0.size(); i++){
                x_next[i] = x[i] + tau * discrepancy(A,b,x);
            }
            x = x_next;
        }
        return x;
    }

    template<std::floating_point T>
    std::vector<T> Gauss_Seidel (const CSR_matrix_space::CSR_matrix<T>& A, const std::vector<T>& b, T discrepancy0, const std::vector<T>& x0){

    }

    template<std::floating_point T>
    friend std::vector<T> Jacobi (const CSR_matrix_space::CSR_matrix<T>& A, const std::vector<T>& b, T discrepancy0, const std::vector<T>& x0){
        std::vector<T> x = x0;
        std::vector<T> x_next = x0;
        while(x_next == x0 || discrepancy(A, b, x) >= discrepancy0){
            for(std::size_t i = 0; i < x0.size(); i++){
                std::vector<T> temp(x0.size());
                for(std::size_t j = 0; j < temp.size(); j++){
                    std::size_t f1 =
                }
            }
        }
    }*/
}

#endif //SLAE_SLAE_SOLVER_HPP
