#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "../src/slae_solver.hpp"

TEST(tridiagonal_tests, tridiagonal_test)
{
    int N = 10;
    Tridiagonal_matrix<double> temp;
    std::vector<double>a(N);
    std::vector<double>b(N);
    std::vector<double>c(N);
    std::vector<double>d(N);
    std::vector<double>x(N);
    std::vector<double>x1;

    std::ifstream fileA;
    std::ifstream filed;
    std::ifstream filex;
    fileA.open("/home/ilya/slae_lab/py_tests/tridiagonal_test/testA.txt");
    filed.open("/home/ilya/slae_lab/py_tests/tridiagonal_test/testd.txt");
    filex.open("/home/ilya/slae_lab/py_tests/tridiagonal_test/testx.txt");
//    filed.open("../../py_tests/tridiagonal_test/testd.txt");
//    filex.open("../../py_tests/tridiagonal_test/testx.txt");
    for (int j = 0; j < 10; j++){
        for(int i = 0; i < N; i++){
            fileA>>a[i];
        }
        for(int i = 0; i < N; i++){
            fileA>>b[i];
        }
        for(int i = 0; i < N; i++){
            fileA>>c[i];
        }
        for(int i = 0; i < N; i++){
            filed>>d[i];
        }
        for(int i = 0; i < N; i++){
            filex>>x[i];
        }
        temp.change_matrix(a,b,c);
        x1 = tridiagonal_matrix_algorithm(temp, d);
        for(std::size_t i = 0; i<x1.size(); i++){
            x1[i] = round(x1[i]*1000)/1000;
            ASSERT_DOUBLE_EQ(x1[i], x[i])<<x1[i]<<" "<<x[i]<<" "<<j<<" "<<i<<"\n";
        }
    }
    fileA.close();
    filed.close();
    filex.close();
}

TEST(constructor_test,copy_constructor_test){
    std::vector<double>a(10);
    std::vector<double>b(10);
    std::vector<double>c(10);
    Tridiagonal_matrix<double> A(a,b,c);
    Tridiagonal_matrix<double> B = A;
    std::cout<<B.size();
}