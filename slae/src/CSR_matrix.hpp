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
        discrepancy(const std::vector<T> &b, const std::vector<T> &x) const {
            std::vector<T> ax = (*this) * x;
            for (std::size_t i = 0; i < x.size(); ++i) {
                ax[i] = b[i] - ax[i];
            }
            return ax;
        }

        T euclid_norm(const std::vector<T> &vec) const {
            return std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), T(0)));
        }

        T scalar_multiplication(const std::vector<T> &vector_1, const std::vector<T> &vector_2) const {
            return std::inner_product(vector_1.begin(), vector_1.end(), vector_2.begin(), T(0));
        }

        std::vector<std::size_t> tau_distribution(std::size_t n) const {
            std::vector<std::size_t> tau_distribution_v(n, 0);
            for (std::size_t i = 2; i <= n; i *= 2) {
                std::size_t j = 0;
                while (j < n) {
                    tau_distribution_v[j + n / i] = i - 1 - tau_distribution_v[j];
                    j += (n * 2) / i;
                }
            }
            return tau_distribution_v;
        }

        std::vector<T> chebyshev_polynomials_solutions(long n) const {
            T cos_b = std::cos(M_PI / (2 * static_cast<T>(n)));
            T sin_b = std::sqrt(1 - std::pow(cos_b, 2));
            //T sin_b = std::sin(M_PI / (2 * static_cast<T>(n)));
            T sin_a = 2 * sin_b * cos_b;
            T cos_a = std::pow(cos_b, 2) - std::pow(sin_b, 2);
            std::vector<T> solutions(n);
            solutions[0] = cos_b;
            for (long i = 0; i + 1 < n; i++) {
                solutions[i + 1] = solutions[i] * cos_a - sin_b * sin_a;
                sin_b = solutions[i] * sin_a + sin_b * cos_a;
            }
            return solutions;
        }
        std::vector<T> one_step_Simmetrical_Gauss_Seidel(const std::vector<T> &b, const std::vector<T> &y) const{
            std::vector<T> temp = y;
            T diag_element;
            for (std::size_t j = 0; j < y.size(); j++) {
                temp[j] = b[j];
                for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                    if (j == col_ind[p]){
                        diag_element = data[p];
                        continue;
                    }
                    temp[j] -= (data[p] * temp[col_ind[p]]);
                }
                temp[j] /= diag_element;
            }

            for (int32_t j = y.size() - 1; j >=0 ; j--) {
                temp[j] = b[j];
                for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                    if (j == static_cast<int32_t>(col_ind[p])){
                        diag_element = data[p];
                        continue;
                    }
                    temp[j] -= (data[p] * temp[col_ind[p]]);
                }
                temp[j] /= diag_element;
            }
            return temp;
        }
        std::vector<T> one_step_SSOR(const std::vector<T> &b, const std::vector<T> &y, T omega) const{
            std::vector<T> temp = y;
            T diag_element;
            for (std::size_t j = 0; j < y.size(); j++) {
                T temp_value = temp[j];
                temp[j] = b[j];
                for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                    if (j == col_ind[p]){
                        diag_element = data[p];
                        continue;
                    }
                    temp[j] -= (data[p] * temp[col_ind[p]]);
                }
                temp[j] *=omega;
                temp[j] /= diag_element;
                temp[j] -= (omega - 1) * temp_value;
            }

            for (int32_t j = y.size() - 1; j >=0 ; j--) {
                T temp_value = temp[j];
                temp[j] = b[j];
                for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                    if (j == static_cast<int32_t>(col_ind[p])){
                        diag_element = data[p];
                        continue;
                    }
                    temp[j] -= (data[p] * temp[col_ind[p]]);
                }
                temp[j] *=omega;
                temp[j] /= diag_element;
                temp[j] -= (omega - 1) * temp_value;
            }
            return temp;
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
                } else if (p + 1 == N) {
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
        std::vector<T> Jacobi(const std::vector<T> &b, const std::vector<T> &x0, T tolerance0) {
            std::vector<T> x = x0;
            T diag_element;
            while (euclid_norm(discrepancy(b, x)) >= tolerance0) {
                std::vector<T> temp(x0.size());
                for (std::size_t j = 0; j < x0.size(); j++) {
                    for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                        if (j == col_ind[p]) {
                            diag_element = data[p];
                            continue;
                        }
                        temp[j] += (data[p] * x[col_ind[p]]);
                    }
                    temp[j] = (b[j] - temp[j])/ diag_element;
                }
                x = temp;
            }
            return x;
        }

        std::vector<T> Gauss_Seidel(const std::vector<T> &b, T tolerance0,
                                    const std::vector<T> &x0) const {
            std::vector<T> x = x0;
            T diag_element;
            while (euclid_norm(discrepancy(b, x)) >= tolerance0) {
                for (std::size_t j = 0; j < x.size(); j++) {
                    x[j] = b[j];
                    for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                        if (j == col_ind[p]){
                            diag_element = data[p];
                            continue;
                        }
                        x[j] -= (data[p] * x[col_ind[p]]);
                    }
                    x[j] /= diag_element;
                }
            }
            return x;
        }

        T estimate_lambda_max(const std::vector<T> &r0, T tolerance0) const {
            std::vector<T> r = r0;
            T lambda, next_lambda;
            T r_err = tolerance0;
            std::vector<T> next_r;

            next_r = (*this) * r;
            T norma = euclid_norm(next_r);
            lambda = scalar_multiplication(next_r, r);
            std::transform(next_r.begin(), next_r.end(), next_r.begin(), [norma](T a) { return a / norma; });
            r = next_r;
            while (r_err >= tolerance0) {
                next_r = (*this) * r;
                norma = euclid_norm(next_r);
                next_lambda = scalar_multiplication(next_r, r);
                std::transform(next_r.begin(), next_r.end(), next_r.begin(), [norma](T a) { return a / norma; });
                r_err = std::abs((lambda - next_lambda)) / lambda;
                lambda = next_lambda;
                r = next_r;
            }
            return lambda;
        }

        std::vector<T> MPI(const std::vector<T> &b, T tau, T tolerance0,
                           const std::vector<T> &x0) const {
            std::vector<T> x = x0;
            std::vector<T> new_discrepancy = discrepancy(b, x);
            while (euclid_norm(new_discrepancy) >= tolerance0) {

                for (std::size_t i = 0; i < x0.size(); i++) {
                    x[i] = x[i] + tau * new_discrepancy[i];
                }
                new_discrepancy = discrepancy(b, x);
            }
            return x;
        }

        std::vector<T>
        chebyshev_MPI(const std::vector<T> &b, T tolerance0, const std::vector<T> &x0, std::size_t n, T lambda_max,
                  T lambda_min) const {
            std::vector<T> x = x0;
            std::vector<std::size_t> tau_distribution_v = tau_distribution(n);
            std::vector<T> tau_v = chebyshev_polynomials_solutions(n);
            std::transform(tau_v.begin(), tau_v.end(), tau_v.begin(), [lambda_max, lambda_min](T a) {
                return 1 / (a * (lambda_max - lambda_min) / 2 + (lambda_max + lambda_min) / 2);
            });
            std::vector<T> new_discrepancy = discrepancy(b, x);
            std::size_t count = 0;
            while (euclid_norm(new_discrepancy) >= tolerance0) {
                if (count >= n) count = 0;
                for (std::size_t i = 0; i < x0.size(); i++) {
                    x[i] = x[i] + tau_v[count] * new_discrepancy[i];
                }
                count++;
                new_discrepancy = discrepancy(b, x);
            }
            return x;
        }
        /*
        std::vector<T>
        quick_MPI_debug(const std::vector<T> &b, T tolerance0, const std::vector<T> &x0, std::size_t n, T lambda_max,
                  T lambda_min) const {
            std::vector<T> x = x0;
            std::vector<std::size_t> tau_distribution_v = tau_distribution(n);
            std::vector<T> tau_v = chebyshev_polynomials_solutions(n);
            std::transform(tau_v.begin(), tau_v.end(), tau_v.begin(), [lambda_max, lambda_min](T a) {
                return 1 / (a * (lambda_max - lambda_min) / 2 + (lambda_max + lambda_min) / 2);
            });
            std::vector<T> new_discrepancy = discrepancy(b, x);
            std::size_t count = 0;
            std::size_t counter = 0;
            while (euclid_norm(new_discrepancy) >= tolerance0) {
                if (abs(x[0]) > 1000) std::cout<<x[0]<<" "<<counter<<" "<<count<<" False \n";
                if (abs(x[1]) > 1000) std::cout<<x[1]<<" "<<counter<<" "<<count<<" False \n";
                if (abs(x[2]) > 1000) std::cout<<x[2]<<" "<<counter<<" "<<count<<" False \n";
                if (count >= n) count = 0;
                for (std::size_t i = 0; i < x0.size(); i++) {
                    x[i] = x[i] + tau_v[count] * new_discrepancy[i];
                }

                count++;
                counter++;
                new_discrepancy = discrepancy(b, x);
            }
            std::cout<<'\n'<<counter<<"\n";
            return x;
        }
         */
        std::vector<T> SOR(const std::vector<T> &b, T tolerance0,
                                    const std::vector<T> &x0, T omega) const {
            std::vector<T> x = x0;
            T diag_element;
            while (euclid_norm(discrepancy(b, x)) >= tolerance0) {
                for (std::size_t j = 0; j < x.size(); j++) {
                    T temp = x[j];
                    x[j] = b[j] ;
                    for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                        if (j == col_ind[p]){
                            diag_element = data[p];
                            continue;
                        }
                        x[j] -= (data[p] * x[col_ind[p]]);
                    }
                    x[j] *=omega;
                    x[j] /= diag_element;
                    x[j] -= (omega - 1) * temp;
                }
            }
            return x;
        }
        std::vector<T> Chebyshev_Simmetrical_Gauss_Seidel(const std::vector<T> &b, T tolerance0,
                           const std::vector<T> &x0, T spectral_radius) const {
            std::vector<T> y_0 = x0;
            std::vector<T> y (x0.size());
            std::vector<T> y_1(x0.size());
            y = one_step_Simmetrical_Gauss_Seidel(b, y_0);


            T mu_0 = T(1);
            T mu = 1/spectral_radius;
            T mu_1;

            while (euclid_norm(discrepancy(b, y_1)) >= tolerance0) {
                mu_1 = 2 * mu / spectral_radius - mu_0;
                y = one_step_Simmetrical_Gauss_Seidel(b, y);
                for(std::size_t i = 0; i < y.size(); i++){
                    y_1[i] = (2 * mu * y[i]/spectral_radius - mu_0 * y_0[i])/mu_1;
                }
                y_0 = y;
                y = y_1;

                mu_0 = mu;
                mu = mu_1;
            }
            return y_1;
        }
        std::vector<T> Chebyshev_SSOR(const std::vector<T> &b, T tolerance0,
                                      const std::vector<T> &x0, T spectral_radius, T omega){
            std::vector<T> y_0 = x0;
            std::vector<T> y (x0.size());
            std::vector<T> y_1(x0.size());
            y = one_step_SSOR(b, y_0, omega);


            T mu_0 = T(1);
            T mu = 1/spectral_radius;
            T mu_1;

            while (euclid_norm(discrepancy(b, y_1)) >= tolerance0) {
                mu_1 = 2 * mu / spectral_radius - mu_0;
                y = one_step_SSOR(b, y, omega);
                for(std::size_t i = 0; i < y.size(); i++){
                    y_1[i] = (2 * mu * y[i]/spectral_radius - mu_0 * y_0[i])/mu_1;
                }
                y_0 = y;
                y = y_1;

                mu_0 = mu;
                mu = mu_1;
            }
            return y_1;
        }
        std::vector<T> steepest_descent(const std::vector<T> &b, T tolerance0, const std::vector<T> &x0){
            std::vector<T> x = x0;
            std::vector<T> discrepancy_v = discrepancy(b,x);
            std::vector<T> ar = (*this) * discrepancy_v;
            T alpha = scalar_multiplication(discrepancy_v, discrepancy_v) / scalar_multiplication(discrepancy_v, ar);
            while (euclid_norm(discrepancy_v) >= tolerance0) {
                x[0] = x[0] + alpha * discrepancy_v[0];
                for (std::size_t i = 1; i < x0.size(); i++) {
                    x[i] = x[i] + alpha * discrepancy_v[i];
                    discrepancy_v[i-1] = discrepancy_v[i - 1] - alpha * ar[i - 1];
                }
                discrepancy_v.back() = discrepancy_v.back() - alpha * ar.back();
                ar = (*this) * discrepancy_v;
                alpha = scalar_multiplication(discrepancy_v, discrepancy_v) / scalar_multiplication(discrepancy_v, ar);
            }
            return x;
        }

        std::vector<T> Polak_balls (const std::vector<T> &b, T tolerance0, const std::vector<T> &x0){
            std::vector<T> x_prev (x0.size());
            std::vector<T> x_next = x0;
            std::vector<T> delta(x0.size());
            std::vector<T> r =  discrepancy(b,x_next);
            std::transform(r.begin(), r.end(), r.begin(), [](auto& a){return -a;});
            std::vector<T> a_r = (*this) * r;
            T r_scalar;
            T r_a_r_scalar;
            T a ;
            T r_delta_scalar;
            std::vector <T> a_delta;
            T delta_a_delta_scalar ;
            T delta_a_r_scalar;
            T r_a_delta_scalar;
            T alpha;
            T betta;
            while (euclid_norm(r) >= tolerance0) {
                std::transform(x_next.begin(), x_next.end(), x_prev.begin(), delta.begin(), std::minus<T>()); //нужно задать вычитание векторов;
                a_r = (*this) * r;
                r_scalar = scalar_multiplication(r, r);
                r_a_r_scalar = scalar_multiplication(r, a_r);
                a = r_scalar/r_a_r_scalar;
                r_delta_scalar = scalar_multiplication(r, delta);
                a_delta = (*this) * delta;
                delta_a_delta_scalar = scalar_multiplication(delta, a_delta);
                delta_a_r_scalar = scalar_multiplication(delta, a_r);
                r_a_delta_scalar = scalar_multiplication(r, a_delta);
//                alpha = (r_scalar * delta_a_delta_scalar - r_delta_scalar * r_a_delta_scalar)/
//                        (r_a_r_scalar*delta_a_delta_scalar - r_a_delta_scalar * r_a_delta_scalar);
//                betta = (alpha * r_a_delta_scalar - r_delta_scalar)/delta_a_delta_scalar;
                betta = (r_scalar * delta_a_r_scalar - r_delta_scalar * r_a_r_scalar)/(delta_a_delta_scalar * r_a_r_scalar - r_a_delta_scalar * r_a_delta_scalar);
                alpha = (betta * r_a_delta_scalar + r_scalar)/r_a_r_scalar;
                *x_prev.begin() = *x_next.begin();
                *x_next.begin() = *x_prev.begin() - alpha * *r.begin() + betta * *delta.begin();
                for (std::size_t i = 1; i < x0.size(); i++) {
                    x_prev[i] = x_next[i];
                    x_next[i] = x_prev[i] - alpha * r[i] + betta * delta[i];
                    r[i-1] = - r[i - 1] + a * a_r[i - 1];
                }
                r.back() = - r.back() + a * a_r.back();
            }
            return x_next;
        }
    };
}
#endif //SLAE_CSR_MATRIX_HPP
