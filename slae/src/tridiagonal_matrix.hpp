//
// Created by ilya on 05.02.23.
//
#ifndef SLAE_TRIDIAGONAL_MATRIX_HPP
#define SLAE_TRIDIAGONAL_MATRIX_HPP

#include <vector>
#include <utility>
#include <algorithm>
#include <concepts>

template <typename T>
concept aritmetical = std::is_floating_point<T>::value || std::is_integral<T>::value;


template <aritmetical T>
class Tridiagonal_matrix{
private:
    template<typename M>
    struct Triads{
        T a;
        T b;
        T c;
        Triads(const T _a,const T _b,const T _c):a(_a), b(_b), c(_c){};
        Triads(const Triads<M> &Tr):a(Tr.a), b(Tr.b), c(Tr.c){};
        Triads() = default;
        Triads& operator = (const Triads<M> &Tr) {
            a = Tr.a;
            b = Tr.b;
            c = Tr.c;
        }
        ~Triads() = default;
    };

    std::vector<Triads<T>> data;
    std::size_t N;
public:
    Tridiagonal_matrix(){
        data.clear();
        N = 0;
    };
    Tridiagonal_matrix(const std::vector<T> &a,const std::vector<T> &b,const std::vector<T> &c){
        data.clear();
        N = a.size();
        for(std::size_t i = 0; i<a.size(); i++){
            data.emplace_back(a[i], b[i],c[i]);
        }
    }
    Tridiagonal_matrix(std::initializer_list<T> A):data(A){};

    Tridiagonal_matrix(const Tridiagonal_matrix<T>& A){
        data.resize(A.size());
        N = A.size();
        for(std::size_t i = 0; i < A.size(); i++){
            data[i] = A[i];
        }
     }
    Tridiagonal_matrix& operator = (const Tridiagonal_matrix<T> &A){
        data.resize(A.size());
        N = A.size();
        for(std::size_t i = 0; i < A.size(); i++){
            data[i] = A[i];
        }
    }
    Tridiagonal_matrix(const Tridiagonal_matrix<T>&& A) noexcept :  Tridiagonal_matrix(std::exchange(A.data, nullptr)){
        N = A.size();
    }

    Tridiagonal_matrix& operator = (Tridiagonal_matrix<T>&& A) noexcept{
        N = A.size();
        std::swap(data, A.data);
    }

    ~Tridiagonal_matrix() = default;
    std::size_t size() const{
        return N;
    }
    void change_matrix(const std::vector<T> &a, const std::vector<T> &b, const std::vector<T> &c)  {
        data.reserve(a.size());
        N = a.size();
        for (std::size_t i = 0; i < a.size(); i++) {
            data[i].a = a[i];
            data[i].b = b[i];
            data[i].c = c[i];
        }
    }
    Triads<T>& operator [](const std::size_t i) {
        return data[i];
    }
    Triads<T> operator [](const std::size_t i) const{
        return data[i];
    }


    T operator()(const std::size_t i,const std::size_t j) const {
        if(j + 1 == i){
            return data[i].a;
        }
        else if(j == i){
            return data[i].b;
        }
        else if(j == i+1){
            return data[i].c;
        }
        throw std::invalid_argument("invalid argument");
    }
    T& operator()(std::size_t i, std::size_t j){
        if(j + 1 == i){
            return data[i].a;
        }
        else if(j == i){
            return data[i].b;
        }
        else if(j == i+1){
            return data[i].c;
        }
        throw std::invalid_argument("invalid argument");
    }

};


#endif //SLAE_TRIDIAGONAL_MATRIX_HPP
