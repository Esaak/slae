//
// Created by ilya on 05.02.23.
//
#ifndef SLAE_TRIDIAGONAL_MATRIX_HPP
#define SLAE_TRIDIAGONAL_MATRIX_HPP

#include <vector>
#include <utility>

namespace {
    template <typename T>
    struct Triads{
    public:
        T a;
        T b;
        T c;
    };
};
template <typename T>
class Tridiagonal_matrix{
private:
    std::vector<Triads<T>> data;
public:
    void create_matrix(const std::vector<T> &a, const std::vector<T> &b, const std::vector<T> &c)  {
        data.clear();
        data.reserve(a.size());
        for (std::size_t i = 0; i < a.size(); i++) {
            data[i].a = a[i];
            data[i].b = b[i];
            data[i].c = c[i];
        }
    }
    Triads<T> operator()(std::size_t i) const {

        return data[i];
    }
    Triads<T>& operator()(std::size_t i){
        return data[i];
    }
};


#endif //SLAE_TRIDIAGONAL_MATRIX_HPP
