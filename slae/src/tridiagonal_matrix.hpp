//
// Created by ilya on 05.02.23.
//
#ifndef SLAE_TRIDIAGONAL_MATRIX_HPP
#define SLAE_TRIDIAGONAL_MATRIX_HPP

#include <vector>
#include <utility>
#include <algorithm>
#include <concepts>
#include <ranges>
template <typename T>
concept aritmetical = std::is_floating_point<T>::value || std::is_integral<T>::value;


template <aritmetical T>
class Tridiagonal_matrix{
private:
    struct Triads{
        T a;
        T b;
        T c;
    };

    std::vector<Triads> data;
    std::size_t N;
public:
    std::size_t size() const{
        return N;
    }
    void change_matrix(const std::vector<T> &a, const std::vector<T> &b, const std::vector<T> &c)  {
        data.reserve(a.size());
        N = a.size();
        //std::ranges::copy(a.begin(), a.end(), );
        //for(auto& it: std::views::zip(a,b, c)){

        for (std::size_t i = 0; i < a.size(); i++) {
            data[i].a = a[i];
            data[i].b = b[i];
            data[i].c = c[i];
        }
    }
    Triads & operator [](const std::size_t i) {
        return data[i];
    }
    const Triads & operator [](const std::size_t i) const{
        return data[i];
    }


    const T& operator()(const std::size_t i,const std::size_t j) const {
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
