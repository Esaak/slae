//
// Created by ilya on 19.04.23.
//
#ifndef SLAE_MATH_HPP
#define SLAE_MATH_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

template<typename T, std::floating_point U>
std::vector<U> discrepancy(const T &A, const std::vector<U> &b, const std::vector<U> &x) {
    std::vector<U> ax = A * x;
    for (std::size_t i = 0; i < x.size(); ++i) {
        ax[i] = b[i] - ax[i];
    }
    return ax;
}
template<std::floating_point T>
T euclid_norm(const std::vector<T> &vec) {
    return std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), T(0)));
}

template<std::floating_point T>
T scalar_multiplication(const std::vector<T> &vector_1, const std::vector<T> &vector_2) {
    return std::inner_product(vector_1.begin(), vector_1.end(), vector_2.begin(), T(0));
}

#endif //SLAE_MATH_HPP
