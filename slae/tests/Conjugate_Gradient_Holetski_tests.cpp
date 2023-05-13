#include <gtest/gtest.h>
#include <string>
#include <fstream>
#include "../src/slae_solver.hpp"


const std::string py_path = "/home/ilya/SLAE/slae/py_tests/Chebyshev_SSOR_tests/";
TEST(Conjugate_Gradient_Holetski_tests, Conjugate_Gradient_Holetski_big_test) {
    std::size_t N = 100;
    std::size_t n = 10;
    std::vector<std::vector<double>> matrix;
    std::vector<double> b(N);
    std::vector<double> x(N);
    std::vector<double> x0(N);
    double tolerance = std::pow(10, -5);
    std::ifstream fileM;
    std::ifstream fileb;
    std::ifstream filex;
    fileM.open(py_path + "test_data.txt");
    fileb.open(py_path + "b.txt");
    filex.open(py_path + "x.txt");
    for (std::size_t k = 0; k < n; k++) {
        matrix.resize(N);
        for (std::size_t i = 0; i < N; i++) {
            std::vector<double> line(N);
            for (std::size_t j = 0; j < N; j++) {
                fileM >> line[j];
            }
            filex >> x[i];
            fileb >> b[i];
            matrix[i] = line;
        }
        Mrx::Matrix<double> M(matrix);
        std::vector<double> result = solvers::Conjugate_Gradient_Holetski(M, b, x0, tolerance);
        for (std::size_t p = 0; p < N; p++) {
            ASSERT_NEAR(result[p], x[p], tolerance);
        }
    }
}