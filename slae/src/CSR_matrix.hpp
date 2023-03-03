//
// Created by ilya on 12.02.23.
//

#ifndef SLAE_CSR_MATRIX_HPP
#define SLAE_CSR_MATRIX_HPP


#include "Matrix.hpp"


namespace DOK_space {
    template<typename T, IsArithmetical<T> = true>
    struct DOK {
        std::size_t i;
        std::size_t j;
        T value;

        bool operator<(const DOK &A) {
            if (i < A.i) {
                return true;
            } else if (i == A.i && j <= A.j) {
                return true;
            }
            return false;
        }
    };
}
namespace CSR_matrix_space {
    template<typename T, IsArithmetical<T> = true>
    class CSR_matrix {
    private:
        std::vector<T> data;
        std::vector<std::size_t> col_ind;
        std::vector<std::size_t> row_indx;
    public:
        void change_matrix(const std::vector<DOK_space::DOK<T>> &A) {
            data.clear();
            std::size_t N = A.size();
            data.resize(N);
            col_ind.resize(N);
            row_indx.resize(A.back().i + 2);
            row_indx[0] = 0;
            std::size_t count = 1;
            for (std::size_t p = 0; p < N; p++) {
                data[p] = A[p].value;
                col_ind[p] = A[p].j;
                if (p > 0 && A[static_cast<int>(p) - 1].i != A[p].i) {
                    std::fill_n(row_indx.begin() + A[static_cast<int>(p) - 1].i + 1,
                                        A[p].i - A[static_cast<int>(p) - 1].i, p);
                    count += A[p].i - A[static_cast<int>(p) - 1].i;
                }
            }
            row_indx.back() = N;
        }

        T &operator()(std::size_t i, std::size_t j) {
            long f1 = static_cast<long>(row_indx[i]);
            long f2 = static_cast<long>(row_indx[i + 1]);
            decltype(col_ind.begin()) result = std::find(col_ind.begin() + f1, col_ind.begin() + f2, j);
            if (result != col_ind.begin() + f2) {
                return data[*result];
            }

            throw std::invalid_argument("invalid argument");
        }

        T operator()(std::size_t i, std::size_t j) const {
            long f1 = static_cast<long>(row_indx[i]);
            long f2 = static_cast<long>(row_indx[i + 1]);
            decltype(col_ind.begin()) result = std::find(col_ind.begin() + f1, col_ind.begin() + f2, j);
            if (result != col_ind.begin() + f2) {
                return data[*result];
            }
            return static_cast<T>(0);
        }

        std::vector<T> operator*(const std::vector<T> &D) {
            std::vector<T> answ(row_indx.size() - 1);
            for (std::size_t i = 0; i < D.size(); i++) {
                if (D[i] == 0) {
                    answ[i] = 0;
                    continue;
                } else {
                    std::size_t f1 = row_indx[i];
                    std::size_t f2 = row_indx[i + 1];
                    for (std::size_t p = f1; p < f2; ++p) {
                        answ[i] += data[p] * D[col_ind[p]];
                    }
                }
            }
            return answ;
        }
    };
}
#endif //SLAE_CSR_MATRIX_HPP
