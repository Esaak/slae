//
// Created by ilya on 09.02.23.
//
#include <iostream>
#include <vector>
#include <gtest/gtest.h>
#include "../tridiagonal_matrix.hpp"

TEST(FirstTest, Test) {
    std::vector<double>a (10);
    std::vector<double>b (10);
    std::vector<double>c (10);
    Tridiagonal_matrix<double>r(a,b,c);
    std::cout<<r(0).a<<" ";
    Tridiagonal_matrix<double> v(r);
    std::cout<<v(0).a;
}