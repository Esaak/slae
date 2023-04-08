#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include "../../src/CSR_matrix.hpp"

using namespace DOK_space;
using namespace CSR_matrix_space;


template<typename T>
void apply_vector(auto &file, std::vector<T> &data) {
    std::string s, temp;
    getline(file, s);
    std::stringstream ss(s);
    while (ss >> temp) {
        std::stringstream str;
        str << temp;
        T f;
        str >> f;
        data.push_back(std::move(f));
    }
}

TEST(TEST_08_04, first_number_MPI_test){
//    double c = 1.0;
//    double b = 30.0;
//    double a = 12.0;
    std::size_t N = 289;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::ifstream fileL_min;
    std::ifstream fileL_max;
    std::vector<double> lambda_max;
    std::vector<double> lambda_min;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b_v;
    std::vector<double> x;
    fileA.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_data.txt");
    filei.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_i.txt");
    filej.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_j.txt");
    fileb.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_b.txt");
    filex.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_x.txt");
    fileL_min.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_L_min.txt");
    fileL_max.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_L_max.txt");
    apply_vector<double>(fileA, data);
    apply_vector<double>(filex, x);
    apply_vector<double>(fileb, b_v);
    apply_vector<std::size_t>(filei, i);
    apply_vector<std::size_t>(filej, j);
    apply_vector<double>(fileL_max, lambda_max);
    apply_vector<double>(fileL_min, lambda_min);
    std::ofstream file;
    std::ofstream file1;
    file.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_count_MPI.txt");
    file1.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_r_MPI.txt");
    std::vector<DOK<double>> D;


    for (std::size_t z = 0; z < data.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    CSR_matrix<double> M(D, N, N);
    double tau = 1 / ((*lambda_max.begin()));
    double tolerance = pow(10, -13);
    std::vector<double> x0(N);
    std::vector<double> result = M.MPI_debug(b_v, tau, tolerance, x0, file, file1);
    for (std::size_t p = 0; p < N; p++) {
        EXPECT_NEAR(result[p], x[p], tolerance);
    }
}
TEST(TEST_08_04, first_number_MPI_opt_test){
//    double c = 1.0;
//    double b = 30.0;
//    double a = 12.0;
    std::size_t N = 289;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::ifstream fileL_min;
    std::ifstream fileL_max;
    std::vector<double> lambda_max;
    std::vector<double> lambda_min;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b_v;
    std::vector<double> x;
    fileA.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_data.txt");
    filei.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_i.txt");
    filej.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_j.txt");
    fileb.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_b.txt");
    filex.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_x.txt");
    fileL_min.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_L_min.txt");
    fileL_max.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_L_max.txt");
    apply_vector<double>(fileA, data);
    apply_vector<double>(filex, x);
    apply_vector<double>(fileb, b_v);
    apply_vector<std::size_t>(filei, i);
    apply_vector<std::size_t>(filej, j);
    apply_vector<double>(fileL_max, lambda_max);
    apply_vector<double>(fileL_min, lambda_min);
    std::ofstream file;
    std::ofstream file1;
    file.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_count_MPI_opt.txt");
    file1.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_r_MPI_opt.txt");
    std::vector<DOK<double>> D;


    for (std::size_t z = 0; z < data.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    CSR_matrix<double> M(D, N, N);
    double tau = 2 / ((*lambda_min.begin()) + (*lambda_max.begin()));
    double tolerance = pow(10, -13);
    std::vector<double> x0(N);
    std::vector<double> result = M.MPI_debug(b_v, tau, tolerance, x0, file, file1);
    for (std::size_t p = 0; p < N; p++) {
        EXPECT_NEAR(result[p], x[p], tolerance);
    }
}
TEST(TEST_08_04, first_number_CHEB_MPI_test){
    std::size_t N = 289;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::ifstream fileL_min;
    std::ifstream fileL_max;
    std::vector<double> lambda_max;
    std::vector<double> lambda_min;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b_v;
    std::vector<double> x;
    fileA.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_data.txt");
    filei.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_i.txt");
    filej.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_j.txt");
    fileb.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_b.txt");
    filex.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_x.txt");
    fileL_min.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_L_min.txt");
    fileL_max.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_L_max.txt");
    apply_vector<double>(fileA, data);
    apply_vector<double>(filex, x);
    apply_vector<double>(fileb, b_v);
    apply_vector<std::size_t>(filei, i);
    apply_vector<std::size_t>(filej, j);
    apply_vector<double>(fileL_max, lambda_max);
    apply_vector<double>(fileL_min, lambda_min);
    std::ofstream file;
    std::ofstream file1;
    file.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_count_CHEB_MPI.txt");
    file1.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_r_CHEB_MPI.txt");
    std::vector<DOK<double>> D;


    for (std::size_t z = 0; z < data.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    CSR_matrix<double> M(D, N, N);
    double tolerance = pow(10, -13);
    std::vector<double> x0(N);
    std::vector<double> result = M.Chebyshev_MPI_debug(b_v, tolerance, x0, 8, lambda_max[0], lambda_min[0], file, file1);
    for (std::size_t p = 0; p < N; p++) {
        EXPECT_NEAR(result[p], x[p], tolerance);
    }
}
TEST(TEST_08_04, first_number_SOR_test){
    std::size_t N = 289;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::ifstream fileL_min;
    std::ifstream fileL_max;
    std::vector<double> lambda_max;
    std::vector<double> lambda_min;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b_v;
    std::vector<double> x;
    fileA.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_data.txt");
    filei.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_i.txt");
    filej.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_j.txt");
    fileb.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_b.txt");
    filex.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_x.txt");
    fileL_min.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_L_min.txt");
    fileL_max.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_L_max.txt");
    apply_vector<double>(fileA, data);
    apply_vector<double>(filex, x);
    apply_vector<double>(fileb, b_v);
    apply_vector<std::size_t>(filei, i);
    apply_vector<std::size_t>(filej, j);
    apply_vector<double>(fileL_max, lambda_max);
    apply_vector<double>(fileL_min, lambda_min);
    std::ofstream file;
    std::ofstream file1;
    file.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_count_SOR.txt");
    file1.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_r_SOR.txt");
    std::vector<DOK<double>> D;


    for (std::size_t z = 0; z < data.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    CSR_matrix<double> M(D, N, N);
    double tolerance = pow(10, -13);
    std::vector<double> x0(N);
    double omega = 1.000345;
    std::vector<double> result = M.SOR_debug(b_v, tolerance, x0, omega, file, file1);
    for (std::size_t p = 0; p < N; p++) {
        EXPECT_NEAR(result[p], x[p], tolerance);
    }
}

TEST(TEST_08_04, first_number_CHEB_MPI_second_test){
    std::size_t N = 289;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::ifstream fileL_min;
    std::ifstream fileL_max;
    std::vector<double> lambda_max;
    std::vector<double> lambda_min;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b_v;
    std::vector<double> x;
    fileA.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_data.txt");
    filei.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_i.txt");
    filej.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_j.txt");
    fileb.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_b.txt");
    filex.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_x.txt");
    fileL_min.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_L_min.txt");
    fileL_max.open("/home/ilya/SLAE/slae/py_tests/TEST_08_04/test_L_max.txt");
    apply_vector<double>(fileA, data);
    apply_vector<double>(filex, x);
    apply_vector<double>(fileb, b_v);
    apply_vector<std::size_t>(filei, i);
    apply_vector<std::size_t>(filej, j);
    apply_vector<double>(fileL_max, lambda_max);
    apply_vector<double>(fileL_min, lambda_min);
    std::ofstream file;
    std::ofstream file1;
    file.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_count_CHEB_MPI_second.txt");
    file1.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_r_CHEB_MPI_second.txt");
    std::vector<DOK<double>> D;


    for (std::size_t z = 0; z < data.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    CSR_matrix<double> M(D, N, N);
    double tolerance = pow(10, -13);
    std::vector<double> x0(N);
    int iterations = 5;
    double lambda_min_scalar = lambda_min[0];
    double lambda_max_scalar = lambda_max[0];
    for(int it = 0; it < iterations; it++){
        std::vector<double> result = M.Chebyshev_MPI_debug_second(b_v, tolerance, x0, 8, lambda_max_scalar, lambda_min_scalar, file, file1);
        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }
        lambda_max_scalar+=1;
        lambda_min_scalar-=1;
    }
}
TEST(TEST_08_04, second_number_MPI_test){
    std::size_t N = 4;
    std::vector<double> data{16.0, 18.0, 21.0, 24.0};
    std::vector<double> i{0, 1, 2, 3};
    std::vector<double> j{0, 1, 2, 3};
    std::vector<double> b{2, 2, 2, 2};
    std::vector<double> c{6, 6, 6, 6};
    double lambda_max =24;
    //double lambda_min =16;
    std::ofstream file;
    std::ofstream file1;
    file.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_count_MPI_coord_ex2.txt");
    file1.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_r_MPI_coord_ex2.txt");
    std::vector<DOK<double>> D;
    for (std::size_t z = 0; z < data.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    CSR_matrix<double> M(D, N, N);
    double tolerance = pow(10, -13);
    std::vector<double> x0(N);
    double tau = 0.9 * 2/lambda_max;
    std::vector<double> result = M.MPI_debug_2(b, tau, tolerance, x0, file, file1);
    std::copy(result.begin(), result.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
}

TEST(TEST_08_04, second_number_MPI_opt_test){
    std::size_t N = 4;
    std::vector<double> data{16.0, 18.0, 21.0, 24.0};
    std::vector<double> i{0, 1, 2, 3};
    std::vector<double> j{0, 1, 2, 3};
    std::vector<double> b{2, 2, 2, 2};
    std::vector<double> c{6, 6, 6, 6};
    double lambda_max =24;
    double lambda_min =16;
    std::ofstream file;
    std::ofstream file1;
    file.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_count_MPI_coord_opt_ex2.txt");
    file1.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_r_MPI_coord_opt_ex2.txt");
    std::vector<DOK<double>> D;
    for (std::size_t z = 0; z < data.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    CSR_matrix<double> M(D, N, N);
    double tolerance = pow(10, -13);
    std::vector<double> x0(N);
    double tau = 2/(lambda_max + lambda_min);
    std::vector<double> result = M.MPI_debug_2(b, tau, tolerance, x0, file, file1);
    std::copy(result.begin(), result.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
}
TEST(TEST_08_04, second_number_Steepest_Descent_test){
    std::size_t N = 4;
    std::vector<double> data{16.0, 18.0, 21.0, 24.0};
    std::vector<double> i{0, 1, 2, 3};
    std::vector<double> j{0, 1, 2, 3};
    std::vector<double> b{2, 2, 2, 2};
    std::vector<double> c{6, 6, 6, 6};
    //double lambda_max =24;
    //double lambda_min =16;
    std::ofstream file;
    std::ofstream file1;
    file.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_count_Steepest_Descent_coord_ex2.txt");
    file1.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_r_Steepest_Descent_coord_ex2.txt");
    std::vector<DOK<double>> D;
    for (std::size_t z = 0; z < data.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    CSR_matrix<double> M(D, N, N);
    double tolerance = pow(10, -13);
    std::vector<double> x0(N);
    //double tau = 2/(lambda_max + lambda_min);
    std::vector<double> result = M.Steepest_Descent_debug_2(b, tolerance, x0, file, file1);
    std::copy(result.begin(), result.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
}

TEST(TEST_08_04, second_number_CHEB_MPI_test){
    std::size_t N = 4;
    std::vector<double> data{16.0, 18.0, 21.0, 24.0};
    std::vector<double> i{0, 1, 2, 3};
    std::vector<double> j{0, 1, 2, 3};
    std::vector<double> b{2, 2, 2, 2};
    std::vector<double> c{6, 6, 6, 6};
    double lambda_max =24;
    double lambda_min =16;
    std::ofstream file;
    std::ofstream file1;
    file.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_count_CHEB_MPI_coord_ex2.txt");
    file1.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_r_CHEB_MPI_coord_ex2.txt");
    std::vector<DOK<double>> D;
    for (std::size_t z = 0; z < data.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    CSR_matrix<double> M(D, N, N);
    double tolerance = pow(10, -13);
    std::vector<double> x0(N);
    //double tau = 2/(lambda_max + lambda_min);
    std::vector<double> result = M.Chebyshev_MPI_debug_2(b, tolerance, x0, 8, lambda_max, lambda_min, file, file1);
    std::copy(result.begin(), result.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
}

TEST(TEST_08_04, second_number_Conjugate_Gradient_test){
    std::size_t N = 4;
    std::vector<double> data{16.0, 18.0, 21.0, 24.0};
    std::vector<double> i{0, 1, 2, 3};
    std::vector<double> j{0, 1, 2, 3};
    std::vector<double> b{2, 2, 2, 2};
    std::vector<double> c{6, 6, 6, 6};

    std::ofstream file;
    std::ofstream file1;
    file.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_count_Conjugate_Gradient_coord_ex2.txt");
    file1.open("/home/ilya/SLAE/slae/slae/tests/TEST_08_04/test_r_Conjugate_Gradient_coord_ex2.txt");
    std::vector<DOK<double>> D;
    for (std::size_t z = 0; z < data.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    CSR_matrix<double> M(D, N, N);
    double tolerance = pow(10, -13);
    std::vector<double> x0(N);

    std::vector<double> result = M.Conjugate_Gradient_debug_2(b, tolerance, x0, file, file1);
    std::copy(result.begin(), result.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
}