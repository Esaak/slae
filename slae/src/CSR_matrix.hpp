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


        std::vector<T>discrepancy(const CSR_matrix_space::CSR_matrix<T>& A, const std::vector<T>& b, const std::vector<T>& x){
            //b- Ax;
            std::vector<T> ax = A * x;
            std::vector<T> answ(b.size());
            for(std::size_t i = 0; i < b.size(); i++){
                answ[i]= b[i] - ax[i];
            }
            return answ;
        }

        T vector_norm(const std::vector<T>& vec){
            return sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), T(0)));
        }


    public:
        CSR_matrix() = default;
        explicit CSR_matrix(const std::vector<DOK_space::DOK<T>> &A, std::size_t row, std::size_t column){ // need , std::size_t row, std::size_t column
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
            }
            row_indx.back() = N;
        }
        T operator()(std::size_t i, std::size_t j) const {
            for(std::size_t p = row_indx[i]; p<row_indx[i + 1]; p++){
                if(col_ind[p] == j){
                    return data[p];
                }
            }
            return static_cast<T>(0);
        }

        std::vector<T> operator*(const std::vector<T> &D) const{
            std::vector<T> answ(row_indx.size() - 1);
            for (std::size_t i = 0; i < D.size(); i++) {
                    std::size_t f1 = row_indx[i];
                    std::size_t f2 = row_indx[i + 1];
                    for (std::size_t p = f1; p < f2; ++p) {
                        answ[i] += data[p] * D[col_ind[p]];
                    }
                }
            return answ;
        }


        std::vector<T> MPI_for3 (const CSR_matrix_space::CSR_matrix<T>& A, const std::vector<T>& b, T tau, T discrepancy0, const std::vector<T>& x0, auto& file){
            std::vector<T> x = x0;
            std::vector<T> x_next = x0;
            std::size_t count = 0;
            T discrepancy_new = vector_norm(discrepancy(A,b,x));
            T delta_discrepancy = discrepancy_new;
            while(discrepancy_new >= discrepancy0 ){
                for(std::size_t i = 0; i < x0.size(); i++){
                    x_next[i] = x[i] + tau * discrepancy(A,b,x)[i];
                }
                x = x_next;
                count+=1;
                if(std::abs(discrepancy_new - vector_norm(discrepancy(A,b,x))) > delta_discrepancy && count>300) break;
                else {
                    delta_discrepancy = std::abs(discrepancy_new - vector_norm(discrepancy(A,b,x)));
                    discrepancy_new = vector_norm(discrepancy(A,b,x));
                }
            }

            file<<count<<" ";

            return x;
        }
        std::vector<T> MPI_for4 (const CSR_matrix_space::CSR_matrix<T>& A, const std::vector<T>& b, T tau, std::size_t iteration_numbers, const std::vector<T>& x0, auto& file){
            std::vector<T> x = x0;
            std::vector<T> x_next = x0;
            std::size_t count = 0;
            T discrepancy_new = vector_norm(discrepancy(A,b,x));
            while(count < iteration_numbers){
                for(std::size_t i = 0; i < x0.size(); i++){
                    x_next[i] = x[i] + tau * discrepancy(A,b,x)[i];
                }
                x = x_next;
                count+=1;
                file<<discrepancy_new<<" ";
                discrepancy_new = vector_norm(discrepancy(A,b,x));
            }
            for(auto& it: x) std::cout<<it<<" ";

            return x;
        }


        std::vector<T> Gauss_Seidel_for4 (const CSR_matrix_space::CSR_matrix<T>& A, const std::vector<T>& b, std::size_t iteration_numbers, const std::vector<T>& x0, auto& file){
            std::vector<T> x = x0;
            std::size_t count = 0;

            T discrepancy_new = vector_norm( discrepancy(A, b, x));
            while( count < iteration_numbers){
                for(std::size_t j = 0; j < x.size(); j++){
                    T t = 0;
                    for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                        if(j == col_ind[p]) continue;
                        t+= (data[p] * x[col_ind[p]])/A(j, j);
                    }
                    x[j] = b[j]/A(j,j) - t;
                }
                file<<discrepancy_new<<" ";
                discrepancy_new = vector_norm( discrepancy(A, b, x));
                count+=1;
            }
            for(auto& it: x) std::cout<<it<<" ";
            return x;
        }


        std::vector<T> Jacobi_for4 (const CSR_matrix_space::CSR_matrix<T>& A, const std::vector<T>& b, std::size_t iteration_numbers, const std::vector<T>& x0, auto& file){
            std::vector<T> x = x0;
            T discrepancy_new = vector_norm( discrepancy(A, b, x));
            std::size_t count = 0;
            while(count < iteration_numbers){
                std::vector<T> temp(x0.size());
                for(std::size_t j = 0; j < temp.size(); j++){
                    for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                        if(j == col_ind[p]) continue;
                        temp[j] += ((data[p] * x[col_ind[p]])/A(j,j));
                    }
                }
                for(std::size_t i = 0; i < x0.size(); i++){
                    x[i] = (b[i]/A(i,i) - temp[i]);
                }
                file<<discrepancy_new<<" ";
                discrepancy_new = vector_norm( discrepancy(A, b, x));
                count+=1;
            }
            for(auto& it: x) std::cout<<it<<" ";
            return x;
        }



        std::vector<T> Jacobi (const CSR_matrix_space::CSR_matrix<T>& A, const std::vector<T>& b, T discrepancy0, const std::vector<T>& x0){
            std::vector<T> x = x0;
            T discrepancy_new = vector_norm( discrepancy(A, b, x));
            while(discrepancy_new < discrepancy0){
                std::vector<T> temp(x0.size());
                for(std::size_t j = 0; j < temp.size(); j++){
                    for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                        if(j == col_ind[p]) continue;
                        temp[j] += (data[p] * x[col_ind[p]])/A(j,j);
                    }
                }
                for(std::size_t i = 0; i < x0.size(); i++){
                    x[i] = (b[i]/A(i,i) - temp[i]);
                }
                discrepancy_new = vector_norm( discrepancy(A, b, x));
            }
            return x;
        }

        std::vector<T> Gauss_Seidel (const CSR_matrix_space::CSR_matrix<T>& A, const std::vector<T>& b, T discrepancy0, const std::vector<T>& x0){
            std::vector<T> x = x0;
            T discrepancy_new = vector_norm( discrepancy(A, b, x));
            while( discrepancy_new < discrepancy0){
                for(std::size_t j = 0; j < x.size(); j++){
                    T t = 0;
                    for (std::size_t p = row_indx[j]; p < row_indx[j + 1]; ++p) {
                        if(j == col_ind[p]) continue;
                        t+= (data[p] * x[col_ind[p]])/A(j, j);
                    }
                    x[j] = b[j]/A(j,j) - t;
                }
                discrepancy_new = vector_norm( discrepancy(A, b, x));
            }
            return x;
        }


        std::vector<T> MPI (const CSR_matrix_space::CSR_matrix<T>& A, const std::vector<T>& b, T tau, T discrepancy0, const std::vector<T>& x0){
            std::vector<T> x = x0;
            std::vector<T> x_next = x0;
            T discrepancy_new = vector_norm(discrepancy(A,b,x));
            while(discrepancy_new >= discrepancy0 ){
                for(std::size_t i = 0; i < x0.size(); i++){
                    x_next[i] = x[i] + tau * discrepancy(A,b,x)[i];
                }
                x = x_next;
                discrepancy_new = vector_norm(discrepancy(A,b,x));
            }
            return x;
        }

    };
}
#endif //SLAE_CSR_MATRIX_HPP
