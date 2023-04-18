#include <gtest/gtest.h>
#include <string>
#include <fstream>
#include "../src/slae_solver.hpp"
#include "../src/Matrix.hpp"

const std::string py_path = "/home/ilya/SLAE/slae/py_tests/GMRES/";

TEST(GMRES_tests, back_gauss_test){
    std::vector<std::vector<double>> data(3);
    data[0] = {2,3,4};
    data[1] = {0,7,6};
    data[2] = {0,0,9};
    std::vector<double> b{3, 2, 4};
    std::vector<double> x0 = {95/126., -2/21., 4/9.};
    Mrx::Matrix<double> M(data);
    std::vector<double> x = solvers::back_Gauss(M,b);
    std::copy(x.begin(), x.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout<<"\n";
    std::copy(x0.begin(), x0.end(), std::ostream_iterator<double>(std::cout, " "));
}

TEST(GMRES_tests, GMRES_test){
    std::vector<std::vector<double>> data(3);
    data[0] = {2,3,4};
    data[1] = {4,7,6};
    data[2] = {5,2,9};
    std::vector<double> b{3, 2, 4};
    std::vector<double> x0 = {-25/8., -1/4., 17/8.};
    Mrx::Matrix<double> M(data);
    std::vector<double> res0(3);
    std::vector<double> x = solvers::GMRES(M,b, res0, 3);
    std::copy(x.begin(), x.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout<<"\n";
    std::copy(x0.begin(), x0.end(), std::ostream_iterator<double>(std::cout, " "));
}

TEST(GMRES_tests, GMRES2_test){
    std::vector<std::vector<double>> data(3);
    data[0] = {-122,301,440};
    data[1] = {4,-72,366};
    data[2] = {510,-324,98};
    std::vector<double> b{3, 2, 4};
    std::vector<double> x0 = {139463/14511703., 67935/14511703., 91139/14511703.};
    Mrx::Matrix<double> M(data);
    std::vector<double> res0(3);
    std::vector<double> x = solvers::GMRES(M,b, res0, 3);
    std::copy(x.begin(), x.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout<<"\n";
    std::copy(x0.begin(), x0.end(), std::ostream_iterator<double>(std::cout, " "));
}

TEST(GMRES_tests, GMRES_big_test){
    std::size_t N = 10;
    std::vector<std::vector<double>> matrix;
    std::vector<double>b(N);
    std::vector<double>x(N);
    std::vector<double>x0(N);
    double tolerance = std::pow(10, -5);
    std::ifstream fileM;
    std::ifstream fileb;
    std::ifstream filex;
    fileM.open(py_path + "matrix.txt");
    fileb.open(py_path + "b.txt");
    filex.open(py_path + "x.txt");
    for(std::size_t k = 0; k < N; k++) {
        matrix.resize(N);
        for (std::size_t i = 0; i < N; i++) {
            std::vector<double> line(N);
            for (std::size_t j = 0; j < N; j++) {
                fileM>>line[j];
            }
            filex>>x[i];
            fileb>>b[i];
            matrix[i] = line;
        }
        Mrx::Matrix<double> M(matrix);
        std::vector<double> result = solvers::GMRES(M, b, x0, N);
        for (std::size_t p = 0; p < N; p++) {
            ASSERT_NEAR(result[p], x[p], tolerance);
        }
        std::cout<<k<<" ";
    }
}