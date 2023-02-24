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
        template<arithmetical T>
        int sgn(T val) {
            return (T(0) < val) - (val < T(0));
        }

        template<arithmetical T>
        std::vector<T> orthogonal_vector(const std::vector<T> &x, std::size_t n) {
            std::vector<T> e(x.size());
            std::vector<T> new_x;
            std::ranges::copy(x.begin(), x.end(), std::back_inserter(new_x));
            e[0] = sgn(new_x[0]) *
                   sqrt(std::inner_product(new_x.begin(), new_x.end(), new_x.begin(), static_cast<T>(0)));
            //e[0] = sqrt(std::inner_product(new_x.begin(), new_x.end(), new_x.begin(), static_cast<T>(0)));
            new_x[0]+=e[0];
            return new_x;
        }

        template<arithmetical T>
        T search_projection(const std::vector<T> &v, const std::vector<T> &a) {
            //на самом деле можно сунуть в один цикл
            T betta = std::inner_product(v.begin(), v.end(), a.begin(), T(0));
            T gamma = std::inner_product(v.begin(), v.end(), v.begin(), T(0));
            return betta / gamma;
        }

        template<arithmetical T>
//здесь неверно, потому что я иду по строке, надо по столбцу
        void multiply_PA(Mrx::Matrix<T> &matrix, const std::vector<T> &v) {
            for (std::size_t j = static_cast<int>(matrix.get_column_size()) - static_cast<int>(v.size()); j < matrix.get_column_size(); j++) {
                T projection = search_projection(v, matrix.get_column(j, matrix.get_column_size() - v.size(), matrix.get_column_size()));
                if (matrix.get_column_size() == j + v.size()) {
                    std::size_t i = j;
                    for (; i +v.size() <= matrix.get_column_size(); i++) {
                        matrix(i, j) -= 2 * v[i-j] * projection;
                    }
                    for (; i < matrix.get_column_size(); i++) {
                        matrix(i, j) = 0;
                    }

                } else {
                    for (std::size_t i = static_cast<int>(matrix.get_column_size()) - static_cast<int>(v.size());
                         i < matrix.get_column_size(); i++) {
                        matrix(i, j) -= 2 * v[i - (static_cast<int>(matrix.get_column_size()) - static_cast<int>(v.size()))] * projection;
                    }
                }
            }
        }
/*
        template<arithmetical T>
        void find_Q2(Mrx::Matrix<T> &Q, const std::vector<T> &v) {
            Mrx::Matrix<T> temp_Q = Q;
            for (std::size_t j = Q.get_column_size() - v.size(); j < Q.get_column_size(); j++) {
                T projection = search_projection(v, temp_Q.get_row(j, Q.get_column_size() - v.size(), Q.get_column_size()));
                //for (std::size_t j = Q.get_column_size() - v.size(); j < Q.get_column_size(); j++) {
                for (std::size_t i = Q.get_column_size() - v.size() ; i < Q.get_column_size(); i++) {
                    Q(i, j) = temp_Q(j, i) - 2 * v[i - (Q.get_column_size() - v.size()) ] * projection;//неверно, нужно смотреть не с нуля

                }
            }
        }
        template<arithmetical T>
        void find_Q3(Mrx::Matrix<T> &Q, const std::vector<T> &v) {
            for (std::size_t i = Q.get_column_size() - v.size(); i < Q.get_column_size(); i++) {
                T projection = search_projection(v, Q.get_row(i, Q.get_column_size() - v.size(), Q.get_column_size()));
                //for (std::size_t j = Q.get_column_size() - v.size(); j < Q.get_column_size(); j++) {
                for (std::size_t j = Q.get_column_size() - v.size() ; j < Q.get_column_size(); j++) {
                    Q(i, j) -= 2 * v[j - (Q.get_column_size() - v.size()) ] * projection;

                }
            }
        }
*/
        template<arithmetical T>
        void find_Q(Mrx::Matrix<T> &Q, const std::vector<T> &v) {
            Mrx::Matrix<T> temp_Q = Q;
            for (std::size_t i = 0; i < Q.get_column_size(); i++) {
                T projection = search_projection(v, temp_Q.get_row(i, Q.get_column_size() - v.size(), Q.get_column_size()));
                //for (std::size_t j = Q.get_column_size() - v.size(); j < Q.get_column_size(); j++) {
                for (std::size_t j = Q.get_column_size() - v.size() ; j < Q.get_column_size(); j++) {
                    Q(i, j) -= 2 * v[j - (Q.get_column_size() - v.size()) ] * projection;//неверно, нужно смотреть не с нуля

                }
            }
        }
    }

    template<arithmetical T>
    std::pair<Mrx::Matrix<T>, Mrx::Matrix<T>> Householder(const Mrx::Matrix<T> &A){
        using namespace _diny;
        std::vector<std::vector<T>> v;
        Mrx::Matrix<T> new_matrix = A;
        Mrx::Matrix<T> Q;

        Q.ones(A.get_column_size(), A.get_row_size());
        for(std::size_t i = 0; i < new_matrix.get_column_size(); ++i){
            std::vector<T> temp1 = new_matrix.get_column(i, i, new_matrix.get_column_size());
            std::vector<T> temp = orthogonal_vector(temp1, i);
            v.emplace_back(temp);
            multiply_PA(new_matrix, v.back());

            find_Q(Q, v.back());
        }

        Q.transponse();
        return std::make_pair(Q, new_matrix);
    }
}

#endif //SLAE_SLAE_SOLVER_HPP
