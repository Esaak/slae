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


        std::vector<T>
        discrepancy(const CSR_matrix_space::CSR_matrix<T> &A, const std::vector<T> &b, const std::vector<T> &x) const {
            std::vector<T> ax = A * x;
            for(std::size_t i = 0; i < x.size(); ++i){
                ax[i]=b[i] - ax[i];
            }
            return ax;
        }

        T euclid_norm(const std::vector<T> &vec) const{
            return std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), T(0)));
        }
        T scalar_multiplication (const std::vector<T> &vector_1, const std::vector<T> &vector_2) const{
            return std::inner_product(vector_1.begin(), vector_1.end(), vector_2.begin(), T(0));
        }

        std::vector<std::size_t> tau_distribution(std::size_t n) const{
            std::vector<std::size_t> tau_distribution_v(n, 0);
            for(std::size_t i = 2; i <= n; i*=2){
                std::size_t j = 0;
                while(j < n){
                    tau_distribution_v[j + n/i] = i - 1 - tau_distribution_v[j];
                    j+=(n * 2)/i;
                }
            }
            return tau_distribution_v;
        }

        std::vector<T> chebyshev_polynomials_solutions(long n) const{
            T cos_b = std::cos(M_PI/(2 * static_cast<T>(n)));
            T sin_b = std::sqrt(1 - std::pow(cos_b,2));
            T sin_a = 2 * sin_b * cos_b;
            T cos_a = std::pow(cos_b,2) - std::pow(sin_b,2);
            std::vector<T> solutions(n);
            solutions[0] = cos_b;
            for(long i = 0; i + 1 < n; i++){
                solutions[i+1] = solutions[i]*cos_a - sin_b * sin_a;
                sin_b = solutions[i] * sin_a + sin_b * cos_a;
            }
            return solutions;
        }

    public:
        CSR_matrix() = default;

        CSR_matrix(const std::vector<DOK_space::DOK<T>> &A, std::size_t row,
                            std::size_t column) { // need , std::size_t row, std::size_t column
            std::size_t N = A.size();
            data.resize(N);
            col_ind.resize(N);
            row_indx.resize(row + 1);
            row_indx[0] = 0;
            for (std::size_t p = 0; p < N; p++) {
                data[p] = A[p].value;
                col_ind[p] = A[p].j;
                if (p > 0 && A[static_cast<int>(p) - 1].i != A[p].i) {
                    std::fill_n(row_indx.begin() + A[static_cast<int>(p) - 1].i + 1,
                                A[p].i - A[static_cast<int>(p) - 1].i, p);
                }
                else if(p + 1==N){
                    int p_indx = std::max(0, static_cast<int>(p) - 1);
                    std::fill(row_indx.begin() + A[p_indx].i + 1, row_indx.end(), N);
                }
            }
            row_indx.back() = N;
        }

        T operator()(std::size_t i, std::size_t j) const {
            for (std::size_t p = row_indx[i]; p < row_indx[i + 1]; p++) {
                if (col_ind[p] == j) {
                    return data[p];
                }
            }
            return static_cast<T>(0);
        }

        std::vector<T> operator*(const std::vector<T> &D) const {
            std::vector<T> answ(row_indx.size() - 1);
            for (std::size_t i = 0; i < D.size(); i++) {
                for (std::size_t p = row_indx[i]; p < row_indx[i + 1]; ++p) {
                    answ[i] += data[p] * D[col_ind[p]];
                }
            }
            return answ;
        }
//        CSR_matrix operator - (const CSR_matrix<T> &other) const{
//            CSR_matrix new_matrix;
//            if()
//
//        }
        std::vector<T> Jacobi(const std::vector<T>&b, const std::vector<T>& x0, T tolerance0) {
            CSR_matrix_space::CSR_matrix<T> A = *this;
            std::vector<T>x = x0;
            std::vector<T> diag_elements;
            std::vector<T> norm_b;
            diag_elements.reserve(x0.size());
            norm_b.reserve(x0.size());

            while (euclid_norm(discrepancy(A, b, x)) >= tolerance0) {
                std::vector<T> temp(x0.size());
                for (std::size_t j = 0; j < x0.size(); j++) {
                    if(diag_elements.size() != x0.size()){
                        diag_elements.emplace_back(A(j,j));
                        norm_b.emplace_back(b[j]/diag_elements[j]);
                    }
                    for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                        if (j == col_ind[p]) {
                            continue;
                        }
                        temp[j] += (data[p] * x[col_ind[p]]) / diag_elements[j];
                    }
                    temp[j] = (norm_b[j] - temp[j]);
                }
                x=temp;
            }
            return x;
        }

        std::vector<T> Gauss_Seidel(const std::vector<T> &b, T tolerance0,
                                    const std::vector<T> &x0) const{
            CSR_matrix_space::CSR_matrix<T> A = *this;
            std::vector<T> x = x0;
            std::vector<T> diag_elements;
            std::vector<T> norm_b;
            diag_elements.reserve(x0.size());
            norm_b.reserve(x0.size());
            while (euclid_norm(discrepancy(A, b, x)) >= tolerance0) {
                for (std::size_t j = 0; j < x.size(); j++) {
                    if(diag_elements.size()!= x0.size()){
                        diag_elements.emplace_back(A(j,j));
                        norm_b.emplace_back(b[j]/diag_elements[j]);
                    }
                    T t = 0;
                    for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                        if (j == col_ind[p]) continue;
                        t += (data[p] * x[col_ind[p]]) / diag_elements[j];
                    }
                    x[j] = norm_b[j] - t;
                }
            }
            return x;
        }

        T estimate_lambda_max (const std::vector<T> &r0, T tolerance0)const{
            CSR_matrix_space::CSR_matrix<T> A = *this;
            std::vector<T> r = r0;
            T lambda, next_lambda;
            T r_err = tolerance0;
            std::vector <T> next_r;

            next_r = A * r;
            T norma = euclid_norm(next_r);
            lambda = scalar_multiplication(next_r, r);
            std::transform(next_r.begin(), next_r.end(), next_r.begin(), [norma](T a){return a/norma;});
            r = next_r;
            while(r_err >= tolerance0){
                next_r = A * r;
                norma = euclid_norm(next_r);
                next_lambda = scalar_multiplication(next_r, r);
                std::transform(next_r.begin(), next_r.end(), next_r.begin(), [norma](T a){return a/norma;});
                r_err = std::abs((lambda - next_lambda))/lambda;
                lambda = next_lambda;
                r = next_r;
            }
            return lambda;
        }

        std::vector<T> MPI(const std::vector<T> &b, T tau, T tolerance0,
                           const std::vector<T> &x0) const{
            CSR_matrix_space::CSR_matrix<T> A = *this;
            std::vector<T> x = x0;
            std::vector<T>new_discrepancy = discrepancy(A, b, x);
            while (euclid_norm(new_discrepancy) >= tolerance0) {

                for (std::size_t i = 0; i < x0.size(); i++) {
                    x[i] = x[i] + tau * new_discrepancy[i];
                }
                new_discrepancy = discrepancy(A, b, x);
            }
            return x;
        }

        std::vector<T> quick_MPI(const std::vector<T>&b, T tolerance0, const std::vector<T>& x0, std::size_t n, T lambda_max, T lambda_min) const{
            std::vector<T> x = x0;
            CSR_matrix A = *this;
            std::vector<std::size_t> tau_distribution_v = tau_distribution(n);
            std::vector<T> tau_v = chebyshev_polynomials_solutions(n);
            std::transform(tau_v.begin(), tau_v.end(), tau_v.begin(), [lambda_max, lambda_min](T a){
                return 1 / (a * (lambda_max - lambda_min)/2 + (lambda_max + lambda_min)/2);
            });
            std::vector<T>new_discrepancy = discrepancy(A, b, x);
            std::size_t count = 0;
            while(euclid_norm(new_discrepancy) >= tolerance0){
                if(count>=n) count = 0;
                for (std::size_t i = 0; i < x0.size(); i++) {
                    x[i] = x[i] + tau_v[count] * new_discrepancy[i];
                }
                count ++;
                new_discrepancy = discrepancy(A, b, x);
            }
            return x;
        }

    };
}
#endif //SLAE_CSR_MATRIX_HPP
