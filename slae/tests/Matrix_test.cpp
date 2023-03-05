//
// Created by ilya on 19.02.23.
//

#include <vector>
#include <gtest/gtest.h>
#include <ranges>
#include <numeric>
#include "../src/Matrix.hpp"
TEST(Matrix_tests, operators_test){
    std::vector<std::vector<int>> test_matrix(3);
    std::vector<std::vector<int>> test_matrix2(3);
    for(std::size_t i = 0; i<3; i++){
        test_matrix[i].resize(3);
        test_matrix2[i].resize(3);
        std::iota(test_matrix[i].begin(), test_matrix[i].end(), i);
    }
    test_matrix2[0][0] = 2;
    test_matrix2[1][1] = 2;
    test_matrix2[2][2] = 2;
    test_matrix[2][2] = 5;
    Mrx::Matrix<int> A(test_matrix);
    Mrx::Matrix<int> B(test_matrix2);
    A*=B;
    for(std::size_t i = 0; i<3;i++){
        for(std::size_t j = 0; j<3;j++){
            std::cout<<A(i, j)<<" ";
        }
        std::cout<<"\n";
    }
}