#include <gtest/gtest.h>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include "../../src/CSR_matrix.hpp"

using namespace DOK_space;
using namespace CSR_matrix_space;

TEST(tests, first){
    std::vector<DOK<double>> D;
    std::vector<int> i{0, 0, 1, 1, 2, 2};
    std::vector<int> j{0, 1, 0, 1, 1, 2};
    std::vector<double> data{10, 1, 1, 7, 0.1, 1};
    std::vector<double> x0{0, 0, 0};
    std::vector<double> b{20, 30, 1};
    std::ofstream file;
    file.open("/home/ilya/SLAE/slae/tests/TEST_04_03/test3_data1.txt");
    for(std::size_t z = 0; z<i.size(); z++){
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    CSR_matrix<double> ANSW(D, 3, 3);
    double r = pow(10, -12);
    double tau = 0.01;
    for(std::size_t p = 0; p < 30; p++) {
        std::vector<double> a = ANSW.MPI_for3(b, tau, r, x0, file);

        for(auto& it: a){
            std::cout<< it<<" ";
        }
        std::cout<<"\n";
        tau += 0.01;
    }
    file.close();
}

TEST(tests, second){
    std::vector<DOK<double>> D;
    std::vector<int> i{0, 0, 0, 1, 1, 1, 2, 2, 2};
    std::vector<int> j{0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> data{12, 17, 3, 17, 15825, 28, 3, 28, 238};
    std::vector<double> x0{1, 1, 1};
    std::vector<double> b{1, 2, 3};
    double c1 = pow(10, 10);
    double c2 = pow(10, -12);
    double c3 = pow(10, 5);
    std::vector<double> C = {c1, c2, c3};
    std::ofstream file_J;
    std::ofstream file_M;
    std::ofstream file_G;
    file_J.open("/home/ilya/SLAE/slae/tests/TEST_04_03/test4_data_jacobi.txt");
    file_M.open("/home/ilya/SLAE/slae/tests/TEST_04_03/test4_data_MPI.txt");
    file_G.open("/home/ilya/SLAE/slae/tests/TEST_04_03/test4_data_Gauss.txt");
    double tau = 0.0001;
    std::size_t iteration_numbers = 10000;
    for(std::size_t z = 0; z<i.size(); z++){
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }

    CSR_matrix<double> ANSW(D, 3, 3);
    ANSW.MPI_for4(b, tau, iteration_numbers, x0, file_M);
    ANSW.Jacobi_for4( b, iteration_numbers, x0, file_J);
    ANSW.Gauss_Seidel_for4(b, iteration_numbers, x0, file_G);
    file_M.close();
    file_J.close();
    file_G.close();
}