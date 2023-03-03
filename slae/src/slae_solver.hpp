//
// Created by ilya on 05.02.23.
//
#ifndef SLAE_SLAE_SOLVER_HPP
#define SLAE_SLAE_SOLVER_HPP

#include "tridiagonal_matrix.hpp"

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
        std::vector<T> orthogonal_vector(const std::vector<T> &x, std::size_t n) {
            T e;
            std::vector<T> new_x;
            std::copy(x.begin(), x.end(), std::back_inserter(new_x));
            e = sgn(new_x[0]) * sqrt(std::inner_product(new_x.begin(), new_x.end(), new_x.begin(), T(0)));
            new_x[0] += e;
            return new_x;
        }

        template<std::floating_point T>
        T search_projection(const std::vector<T> &v, const std::vector<T> &a, const T norm) {
            T betta = std::inner_product(v.begin(), v.end(), a.begin(), T(0));
            return betta / norm;
        }

        template<std::floating_point T>
        void find_R(Mrx::Matrix<T> &matrix, const std::vector<T> &v) {
            std::size_t n = static_cast<int>(matrix.get_column_size()) - static_cast<int>(v.size());
            T norm = std::inner_product(v.begin(), v.end(), v.begin(), T(0));
            for (std::size_t j = n; j < matrix.get_column_size(); j++) {
                T projection = search_projection(v, matrix.get_column(j, n, matrix.get_column_size()), norm);
                if (n == j) {
                    std::size_t i = j;
                    for (; i <= n; i++) {
                        matrix(i, j) -= 2 * v[i - j] * projection;
                    }
                    for (; i < matrix.get_column_size(); i++) {
                        matrix(i, j) = 0;
                    }

                } else {
                    for (std::size_t i = n; i < matrix.get_column_size(); i++) {
                        matrix(i, j) -= 2 * v[i - n] * projection;
                    }
                }
            }
        }

/*

*/
        template<std::floating_point T>
        void find_Q(Mrx::Matrix<T> &Q, const std::vector<T> &v) {
            Mrx::Matrix<T> temp_Q = Q;
            T norm = std::inner_product(v.begin(), v.end(), v.begin(), T(0));
            std::size_t n = static_cast<int>(Q.get_column_size()) - static_cast<int>(v.size());
            for (std::size_t i = 0; i < Q.get_column_size(); i++) {
                T projection = search_projection(v, temp_Q.get_column(i, n, Q.get_column_size()), norm);
                for (std::size_t j = n; j < Q.get_column_size(); j++) {
                    Q(j, i) -= 2 * v[j - n] * projection;
                }
            }
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
        temp.eye(N, N);
        for (std::size_t p = 0; p < matrix.size() - 1; p++) {
            Mrx::Matrix<T> temp1;
            Mrx::Matrix<T> temp2;
            temp1.eye(N, N);
            temp2.eye(matrix[p].size(), matrix[p].size());
            std::size_t i;
            T norm = std::inner_product(matrix[p].begin(), matrix[p].end(), matrix[p].begin(), T(0));
            for (i = N - matrix[p].size(); i < N; ++i) {
                T projection = search_projection(matrix[p], temp2.get_column(i - (N - matrix[p].size()), 0,
                                                                             temp2.get_column_size()), norm);
                for (std::size_t j = N - matrix[p].size(); j < N; j++) {
                    temp1(j, i) -=
                            2 * matrix[p][j - (N - matrix[p].size())] * projection;//неверно, нужно смотреть не с нуля
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

        Q.eye(A.get_column_size(), A.get_row_size());
        for (std::size_t i = 0; i + 1 < new_matrix.get_column_size(); ++i) {
            std::vector<T> temp1 = new_matrix.get_column(i, i, new_matrix.get_column_size());
            std::vector<T> temp = orthogonal_vector(temp1, i);
            find_R(new_matrix, temp);
            find_Q(Q, temp);
        }
        Q.transponse();
        return std::make_pair(Q, new_matrix);
    }
}

#endif //SLAE_SLAE_SOLVER_HPP
