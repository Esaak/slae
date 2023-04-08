#include <gtest/gtest.h>
#include "../src/slae_solver.hpp"
#include <numeric>
#include <utility>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cmath>

TEST(QR_tests, simple_test){
    //int N = 9;
    //int vec_size = 9;
    std::vector<std::vector<double>> test{{1,2,3}, {4,5,6}, {7,8,9}};
    Mrx::Matrix<double> matrix(test);

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
    D.first*=D.second;
    for(std::size_t i = 0; i < 3; i++){
        for(std::size_t j = 0; j < 3; j++){
            std::cout<< D.first(i,j)<<" ";
        }
        std::cout<<"\n";
    }
}

TEST(QR_tests, main_test){
    std::size_t N = 10;
    std::vector<std::vector<double>> matrix;
    std::vector<std::vector<double>> R;
    std::vector<std::vector<double>> Q;
    std::ifstream fileM;
    std::ifstream fileQ;
    std::ifstream fileR;
    fileM.open("/home/ilya/slae_lab/py_tests/QR_test/matrix.txt");
    fileQ.open("/home/ilya/slae_lab/py_tests/QR_test/Q.txt");
    fileR.open("/home/ilya/slae_lab/py_tests/QR_test/R.txt");
    for(std::size_t i = 0; i < N; i++){
        for(std::size_t p = 0; p < N; p++){
            std::vector<double> tm(N);
            std::vector<double> tR(N);
            std::vector<double> tQ(N);
            for(std::size_t q = 0; q < N; q++){
                fileM>>tm[q];
                fileQ>>tQ[q];
                fileR>>tR[q];
            }
            matrix.push_back(tm);
            Q.push_back(tQ);
            R.push_back(tR);
        }
        Mrx::Matrix<double>R_matrix(matrix);
        Mrx::Matrix<double>R_Q(Q);
        Mrx::Matrix<double>R_R(R);
        std::pair<Mrx::Matrix<double>, Mrx::Matrix<double>> answ;
        answ = solvers::Householder(R_matrix);
        for(std::size_t x = 0; x < N; x++){
            for (std::size_t y = 0; y < N; ++y) {
                EXPECT_NEAR(answ.first[x][y], Q[x][y], 0.001)<<x<<" "<<y<<" "<<i<<"\n";

            }
        }

        for(std::size_t x = 0; x < N; x++){
            for (std::size_t y = 0; y < N; ++y) {
                EXPECT_NEAR(answ.second[x][y], R[x][y], 0.001);
            }
        }
        matrix.clear();
        Q.clear();
        R.clear();
    }
}