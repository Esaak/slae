#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "../src/CSR_matrix.hpp"

TEST(CRS_matrix_test, matrix_filling_test){
    std::size_t N = 10;
    std::vector<double>data;
    std::vector<std::size_t>indices;
    std::vector<std::size_t>indptr;
    std::vector<std::size_t>i;
    std::vector<std::size_t>j;
    std::ifstream fileA;
    std::ifstream filed;
    std::ifstream filex;
    std::ifstream filei;
    std::ifstream filej;
    fileA.open("/home/ilya/slae_lab/py_tests/SRC_tests/test_data.txt");
    filed.open("/home/ilya/slae_lab/py_tests/SRC_tests/test_indices.txt");
    filex.open("/home/ilya/slae_lab/py_tests/SRC_tests/test_indptr.txt");
    filei.open("/home/ilya/slae_lab/py_tests/SRC_tests/test_i.txt");
    filej.open("/home/ilya/slae_lab/py_tests/tridiagonal_test/test_j.txt");
    std::string str;
    data = {8.943, 4.212, 6.965, 3.956, 6.448, 8.476, 1.279, 3.074, 2.49, 2.31, 3.842, 1.086, 3.214, 1.574, 3.556, 3.948 };
    indices = {3, 6, 0, 8, 9, 6, 0, 0, 2, 8, 2, 3, 0, 6, 4, 9};
    indptr = {0, 1, 2, 5, 6, 7, 10, 12, 14, 14, 16};
    i = {0, 1, 2, 2, 2, 3, 4, 5, 5, 5, 6, 6, 7, 7, 9, 9};
    j = {3, 6, 0, 8, 9, 6, 0, 0, 2, 8, 2, 3, 0, 6, 4, 9 };
    std::reverse(i.begin(), i.end());
    std::reverse(j.begin(), j.end());
    std::reverse(data.begin(), data.end());
    std::vector<DOK<double>>D;
    for(int z =0; z<data.size(); z++){
        D.emplace_back(DOK<double>{i[z], j[z], data[z]});
    }
    std::sort(D.begin(),D.end());
    CSR_matrix<double> M(D);
    std::reverse(data.begin(), data.end());
    //std::cout<<M.compare(data, indices, indptr);
//    for(std::size_t x = 0; x<i.size(); x++){
//            std::cout<<M(i[x], j[x])<<" ";
//    }
    std::vector<double> d{1,1,1,1,1,1,1,1,1,1};
    std::vector<double> w = M*d;
    std::copy(w.begin(), w.end(), std::ostream_iterator<double>(std::cout, " "));

}