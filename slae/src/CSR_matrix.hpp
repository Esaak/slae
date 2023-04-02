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

        std::vector<T> one_step_Simmetrical_Gauss_Seidel(const std::vector<T> &b, const std::vector<T> &y) const {
            std::vector<T> temp = y;
            T diag_element = (1);
            for (std::size_t j = 0; j < y.size(); j++) {
                temp[j] = b[j];
                for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                    if (j == col_ind[p]) {
                        diag_element = data[p];
                        continue;
                    }
                    temp[j] -= (data[p] * temp[col_ind[p]]);
                }
                temp[j] /= diag_element;
            }

            for (int32_t j = y.size() - 1; j >= 0; j--) {
                temp[j] = b[j];
                for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                    if (j == static_cast<int32_t>(col_ind[p])) {
                        diag_element = data[p];
                        continue;
                    }
                    temp[j] -= (data[p] * temp[col_ind[p]]);
                }
                temp[j] /= diag_element;
            }
            return temp;
        }

        std::vector<T> one_step_SSOR(const std::vector<T> &b, const std::vector<T> &y, T omega) const {
            std::vector<T> temp = y;
            T diag_element = (1);
            for (std::size_t j = 0; j < y.size(); j++) {
                T temp_value = temp[j];
                temp[j] = b[j];
                for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                    if (j == col_ind[p]) {
                        diag_element = data[p];
                        continue;
                    }
                    temp[j] -= (data[p] * temp[col_ind[p]]);
                }
                temp[j] *= omega;
                temp[j] /= diag_element;
                temp[j] -= (omega - 1) * temp_value;
            }

            for (int32_t j = y.size() - 1; j >= 0; j--) {
                T temp_value = temp[j];
                temp[j] = b[j];
                for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                    if (j == static_cast<int32_t>(col_ind[p])) {
                        diag_element = data[p];
                        continue;
                    }
                    temp[j] -= (data[p] * temp[col_ind[p]]);
                }
                temp[j] *= omega;
                temp[j] /= diag_element;
                temp[j] -= (omega - 1) * temp_value;
            }
            return temp;
        }

        std::pair<T, T> math_for_Polak_balls(const std::vector<T> &a_x_prev, const std::vector<T> &a_x_next,
                                             const std::vector<T> &x_prev, const std::vector<T> &x_next,
                                             const std::vector<T> &delta, const std::vector<T> &r) const {
            std::vector<T> a_delta(x_next.size());
            std::transform(a_x_next.begin(), a_x_next.end(), a_x_prev.begin(), a_delta.begin(), std::minus<T>());
            std::vector<T> a_r = (*this) * r;
            T r_scalar = scalar_multiplication(r, r);
            T r_a_r_scalar = scalar_multiplication(r, a_r);
            T r_delta_scalar = scalar_multiplication(r, delta);
            T delta_a_delta_scalar = scalar_multiplication(delta, a_delta);
            T r_a_delta_scalar = scalar_multiplication(r, a_delta);
            T betta = (-r_scalar * r_a_delta_scalar + r_delta_scalar * r_a_r_scalar) /
                      (delta_a_delta_scalar * r_a_r_scalar - r_a_delta_scalar * r_a_delta_scalar);
            T alpha = (-betta * r_a_delta_scalar + r_scalar) / r_a_r_scalar;
            return std::make_pair(alpha, betta);
        }

    public:
        CSR_matrix() = default;

        CSR_matrix(const std::vector<DOK_space::DOK<T>> &A, std::size_t row,
                   std::size_t column) {
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

        std::vector<T> Jacobi(const std::vector<T> &b, const std::vector<T> &x0, T tolerance0) {
            std::vector<T> x = x0;
            T diag_element = T(1);
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
                    temp[j] = (b[j] - temp[j]) / diag_element;
                }
                x = temp;
            }
            return x;
        }

        std::vector<T> Gauss_Seidel(const std::vector<T> &b, T tolerance0,
                                    const std::vector<T> &x0) const {
            std::vector<T> x = x0;
            T diag_element = T(1);
            while (euclid_norm(discrepancy(b, x)) >= tolerance0) {
                for (std::size_t j = 0; j < x.size(); j++) {
                    x[j] = b[j];
                    for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                        if (j == col_ind[p]) {
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

        std::vector<T> Chebyshev_Simmetrical_Gauss_Seidel(const std::vector<T> &b, T tolerance0,
                                                          const std::vector<T> &x0, T spectral_radius) const {
            std::vector<T> y_prev = x0;
            std::vector<T> y_next = one_step_Simmetrical_Gauss_Seidel(b, y_prev);

            T mu_prev = T(1);
            T mu = 1 / spectral_radius;
            T mu_next;

            while (euclid_norm(discrepancy(b, y_next)) >= tolerance0) {
                mu_next = 2 * mu / spectral_radius - mu_prev;
                y_next = one_step_Simmetrical_Gauss_Seidel(b, y_next);
                for (std::size_t i = 0; i < y_next.size(); i++) {
                    y_next[i] = (2 * mu * y_next[i] / spectral_radius - mu_prev * y_prev[i]) / mu_next;
                }
                y_prev = y_next;
                mu_prev = mu;
                mu = mu_next;
            }
            return y_next;
        }

        std::vector<T> MPI(const std::vector<T> &b, T tau, T tolerance0,
                           const std::vector<T> &x0) const {
            std::vector<T> x = x0;
            std::vector<T> r = discrepancy(b, x);
            while (euclid_norm(r) >= tolerance0) {

                for (std::size_t i = 0; i < x0.size(); i++) {
                    x[i] = x[i] + tau * r[i];
                }
                r = discrepancy(b, x);
            }
            return x;
        }

        std::vector<T>
        Chebyshev_MPI(const std::vector<T> &b, T tolerance0, const std::vector<T> &x0, std::size_t n, T lambda_max,
                      T lambda_min) const {
            std::vector<T> x = x0;
            std::size_t count = 0;
            std::vector<std::size_t> tau_distribution_v = tau_distribution(n);
            std::vector<T> tau_v = chebyshev_polynomials_solutions(n);
            std::transform(tau_v.begin(), tau_v.end(), tau_v.begin(), [lambda_max, lambda_min](T a) {
                return 2 / (a * (lambda_max - lambda_min) + (lambda_max + lambda_min));
            });
            std::vector<T> r = discrepancy(b, x);
            while (euclid_norm(r) >= tolerance0) {
                if (count >= n) count = 0;
                for (std::size_t i = 0; i < x0.size(); i++) {
                    x[i] = x[i] + tau_v[count] * r[i];
                }
                count++;
                r = discrepancy(b, x);
            }
            return x;
        }

        std::vector<T> SOR(const std::vector<T> &b, T tolerance0,
                           const std::vector<T> &x0, T omega) const {
            std::vector<T> x = x0;
            T diag_element = (1);
            while (euclid_norm(discrepancy(b, x)) >= tolerance0) {
                for (std::size_t j = 0; j < x.size(); j++) {
                    T temp = x[j];
                    x[j] = b[j];
                    for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                        if (j == col_ind[p]) {
                            diag_element = data[p];
                            continue;
                        }
                        x[j] -= (data[p] * x[col_ind[p]]);
                    }
                    x[j] *= omega;
                    x[j] /= diag_element;
                    x[j] -= (omega - 1) * temp;
                }
            }
            return x;
        }

        std::vector<T> Chebyshev_SSOR(const std::vector<T> &b, T tolerance0,
                                      const std::vector<T> &x0, T spectral_radius, T omega) {
            std::vector<T> y_prev = x0;
            std::vector<T> y_next = one_step_SSOR(b, y_prev, omega);

            T mu_prev = T(1);
            T mu = 1 / spectral_radius;
            T mu_next;

            while (euclid_norm(discrepancy(b, y_next)) >= tolerance0) {
                mu_next = 2 * mu / spectral_radius - mu_prev;
                y_next = one_step_SSOR(b, y_next, omega);
                for (std::size_t i = 0; i < y_next.size(); i++) {
                    y_next[i] = (2 * mu * y_next[i] / spectral_radius - mu_prev * y_prev[i]) / mu_next;
                }
                y_prev = y_next;
                mu_prev = mu;
                mu = mu_next;
            }
            return y_next;
        }

        std::vector<T> Steepest_Descent(const std::vector<T> &b, T tolerance0, const std::vector<T> &x0) {
            std::vector<T> x = x0;
            std::vector<T> r = discrepancy(b, x);
            std::vector<T> a_r = (*this) * r;
            T alpha = scalar_multiplication(r, r) / scalar_multiplication(r, a_r);
            while (euclid_norm(r) >= tolerance0) {
                x[0] = x[0] + alpha * r[0];
                for (std::size_t i = 1; i < x0.size(); i++) {
                    x[i] = x[i] + alpha * r[i];
                    r[i - 1] = r[i - 1] - alpha * a_r[i - 1];
                }
                r.back() = r.back() - alpha * a_r.back();
                a_r = (*this) * r;
                alpha = scalar_multiplication(r, r) / scalar_multiplication(r, a_r);
            }
            return x;
        }

        std::vector<T> Polak_Balls(const std::vector<T> &b, T tolerance0, const std::vector<T> &x0) {
            std::vector<T> x_prev(x0.size()), delta(x0.size());
            std::vector<T> x_next = x0;
            std::vector<T> a_x_prev(x0.size()), a_x_next;
            std::vector<T> r(x0.size());
            std::pair<T, T> iteration_coefficient;
            a_x_next = (*this) * x_next;
            std::transform(b.begin(), b.end(), a_x_next.begin(), r.begin(), std::minus<T>());

            while (euclid_norm(r) >= tolerance0) {
                std::transform(x_next.begin(), x_next.end(), x_prev.begin(), delta.begin(), std::minus<T>());
                iteration_coefficient = math_for_Polak_balls(a_x_prev, a_x_next, x_prev, x_next, delta, r);
                for (std::size_t i = 0; i < x0.size(); i++) {
                    x_prev[i] = x_next[i];
                    x_next[i] =
                            x_prev[i] + iteration_coefficient.first * r[i] + iteration_coefficient.second * delta[i];
                }
                a_x_prev = a_x_next;
                a_x_next = (*this) * x_next;
                std::transform(b.begin(), b.end(), a_x_next.begin(), r.begin(), std::minus<T>());
            }

            return x_next;
        }

        std::vector<T> Conjugate_Gradient(const std::vector<T>& b,T tolerance0, const std::vector<T>& x0) const{
            std::vector<T>r, d, x = x0;
            r = discrepancy(b, x0);
            d = r;
            while(euclid_norm(r)>= tolerance0){
                T r_d_scalar = scalar_multiplication(r, d);
                T alpha = r_d_scalar/ scalar_multiplication(d, (*this) * d);
                for(std::size_t i = 0; i < x.size(); i++){
                    x[i] += alpha * d[i];
                }
                r = discrepancy(b,x);
                T r_scalar = scalar_multiplication(r,r);
                T coef = r_scalar/r_d_scalar;
                std::transform(d.begin(), d.end(), r.begin(), d.begin(), [&coef](auto&& a, auto&&b){
                    return b + coef * a;
                });
            }
            return x;
        }
    };
}
#endif //SLAE_CSR_MATRIX_HPP
