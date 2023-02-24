#include <gtest/gtest.h>
#include "../src/Matrix.hpp"
#include "../src/slae_solver.hpp"
#include<ranges>
#include <numeric>
#include <utility>
#include <cmath>
TEST(QR_tests, main_test){
    //int N = 9;
    //int vec_size = 9;
    std::vector<std::vector<double>> test{{1,2,3}, {4,5,6}, {7,8,9}};
    Mrx::Matrix<double> matrix;
    matrix.change_matrix(test);
    /*
    for(auto&& it: test){
        for(auto&& t:it){
            std::cout<<t<<" ";
        }
        std::cout<<"\n";
    }
    */
    std::pair<Mrx::Matrix<double>, Mrx::Matrix<double>> D;
    D = solvers::Householder(matrix);
    for(std::size_t i = 0; i < 3; i++){
        for(std::size_t j = 0; j < 3; j++){
            std::cout<< D.first(i,j)<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";
    D.second*=D.first;
    for(std::size_t i = 0; i < 3; i++){
        for(std::size_t j = 0; j < 3; j++){
            std::cout<< D.second(i,j)<<" ";
        }
        std::cout<<"\n";
    }
}