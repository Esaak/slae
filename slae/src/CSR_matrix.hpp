//
// Created by ilya on 12.02.23.
//

#ifndef SLAE_CSR_MATRIX_HPP
#define SLAE_CSR_MATRIX_HPP
#include <vector>
#include <concepts>
#include <compare>
#include <cstring>
#include <utility>
#include <algorithm>

template <typename T>
concept aritmetical = std::is_floating_point<T>::value || std::is_integral<T>::value;

template<aritmetical T>
struct DOK{
    std::size_t i;
    std::size_t j;
    T value;
    bool operator < (const DOK& A){
        if(i < A.i){
            return true;
        }
        else if(i<=A.i && j <= A.j){
            return true;
        }
        return false;
    }
};

template<aritmetical T>
class CSR_matrix{
private:
    std::vector<T> data;
    std::vector<std::size_t> col_ind;
    std::vector<std::size_t> row_indx;
public:
    /*
    CSR_matrix() = default;
    ~CSR_matrix() = default;

    CSR_matrix(const std::vector<DOK<T>>& A){

    }
    CSR_matrix(const std::vector<T> &_data, const std::vector<T>& _col, const std::vector<T>& _row){
        data = _data;
        col_ind = _col;
        row_indx = _row;
    }
    CSR_matrix(const CSR_matrix<T>& A){
        data = A.data;
        col_ind = A.col_ind;
        row_indx = A.row_indx;
    }
    //CSR_matrix(std::initializer_list<T> A):data(A){}/тут не определено значение N

    CSR_matrix& operator = (const CSR_matrix<T>& A){
        *this(A);
    }
    CSR_matrix (const CSR_matrix<T>&& A) noexcept {
        *this = std::move(A);
    }
    CSR_matrix& operator = (CSR_matrix<T>&& A) noexcept{
        std::swap(data, A);
    }
    */
    void change_matrix(const std::vector<DOK<T>>& A){
        data.clear();
        std::size_t N = A.size();
        data.resize(N);
        col_ind.resize(N);
        row_indx.resize(A.back().i + 2);
        row_indx[0] = 0;
        std::size_t count = 0;
        for(std::size_t p = 0; p < N; p++){
            data[p] = A[p].value;
            col_ind[p] = A[p].j;
            if(p > 0 && A[static_cast<int>(p) - 1].i != A[p].i){
                for(std::size_t vn = A[static_cast<int>(p) - 1].i; vn < A[p].i; vn++){
                    row_indx[++count] = p;
                }
            }
        }
        row_indx.back() = N;
    }

    T& operator() (std::size_t i, std::size_t j){
        std::size_t f1 = row_indx[i];
        std::size_t f2 = row_indx[i+1];
        for(std::size_t p = f1; p<f2; ++p){
            if(col_ind[p] == j){
                return data[p];
            }
        }
        throw std::invalid_argument("invalid argument");
    }
    T operator() (std::size_t i, std::size_t j) const{
        std::size_t f1 = row_indx[i];
        std::size_t f2 = row_indx[i+1];
        for(std::size_t p = f1; p<f2; ++p){
            if(col_ind[p] == j){
                return data[p];
            }
        }
        return static_cast<T>(0);
    }

    std::vector<T> operator *(const std::vector<T> &D){
        std::vector<T> answ(row_indx.size() - 1);
        for(std::size_t i = 0; i<D.size(); i++){
            if(D[i] == 0){
                answ[i] = 0;
                continue;
            }
            else{
                std::size_t f1 = row_indx[i];
                std::size_t f2 = row_indx[i+1];
                for(std::size_t p = f1; p<f2; ++p){
                    answ[i] += data[p] * D[col_ind[p]];
                }
            }
        }
        return answ;
    }
};
#endif //SLAE_CSR_MATRIX_HPP
