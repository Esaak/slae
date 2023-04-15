#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include "../src/CSR_matrix.hpp"

using namespace DOK_space;
using namespace CSR_matrix_space;

const std::string py_path = "/home/ilya/SLAE/slae/py_tests/";

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

TEST(Iteration_tests, MPI_test) {
    std::size_t N = 100;
    std::size_t n = 10;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::ifstream fileL;
    std::vector<double> lambda;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b;
    std::vector<double> x;
    fileA.open(py_path + "Iteration_tests/test_data.txt");
    filei.open(py_path + "Iteration_tests/test_i.txt");
    filej.open(py_path + "Iteration_tests/test_j.txt");
    fileb.open(py_path + "Iteration_tests/test_b.txt");
    filex.open(py_path + "Iteration_tests/test_x.txt");
    fileL.open(py_path + "Iteration_tests/test_L.txt");

    for (std::size_t it = 0; it < n; it++) {
        apply_vector<double>(fileA, data);
        apply_vector<double>(filex, x);
        apply_vector<double>(fileb, b);
        apply_vector<std::size_t>(filei, i);
        apply_vector<std::size_t>(filej, j);
        apply_vector<double>(fileL, lambda);
        std::vector<DOK<double>> D;
        for (std::size_t z = 0; z < data.size(); z++) {
            D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        }
        CSR_matrix<double> M(D, N, N);
        double tau = 2 / (100 * (*lambda.begin()));
        double tolerance = pow(10, -5);
        std::vector<double> x0(N);

        std::vector<double> result = M.MPI(b, tau, tolerance, x0);

        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        data.clear();
        x.clear();
        b.clear();
        i.clear();
        j.clear();
        lambda.clear();
    }
}

TEST(Iteration_tests, Jacobi_test) {
    std::size_t N = 100;
    std::size_t n = 10;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b;
    std::vector<double> x;
    std::vector<double> result;
    fileA.open(py_path + "Iteration_tests/test_data.txt");
    filei.open(py_path + "Iteration_tests/test_i.txt");
    filej.open(py_path + "Iteration_tests/test_j.txt");
    fileb.open(py_path + "Iteration_tests/test_b.txt");
    filex.open(py_path + "Iteration_tests/test_x.txt");

    for (std::size_t it = 0; it < n; it++) {
        apply_vector<double>(fileA, data);
        apply_vector<double>(filex, x);
        apply_vector<double>(fileb, b);
        apply_vector<std::size_t>(filei, i);
        apply_vector<std::size_t>(filej, j);
        std::vector<DOK<double>> D;
        for (std::size_t z = 0; z < data.size(); z++) {
            D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        }
        CSR_matrix<double> M(D, N, N);
        double tolerance = pow(10, -10);
        data.clear();
        data.resize(N, 0);
        std::vector<double> x0(N);
        result = M.Jacobi(b, tolerance, data);
        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }
        data.clear();
        x.clear();
        b.clear();
        i.clear();
        j.clear();
    }
}

TEST(Iteration_tests, Gauss_Seidel_test) {
    std::size_t N = 100;
    std::size_t n = 10;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b;
    std::vector<double> x;
    fileA.open(py_path + "Iteration_tests/test_data.txt");
    filei.open(py_path + "Iteration_tests/test_i.txt");
    filej.open(py_path + "Iteration_tests/test_j.txt");
    fileb.open(py_path + "Iteration_tests/test_b.txt");
    filex.open(py_path + "Iteration_tests/test_x.txt");

    for (std::size_t it = 0; it < n; it++) {
        apply_vector<double>(fileA, data);
        apply_vector<double>(filex, x);
        apply_vector<double>(fileb, b);
        apply_vector<std::size_t>(filei, i);
        apply_vector<std::size_t>(filej, j);
        std::vector<DOK<double>> D;
        for (std::size_t z = 0; z < data.size(); z++) {
            D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        }
        CSR_matrix<double> M(D, N, N);
        double tolerance = pow(10, -10);
        std::vector<double> x0(N);
        std::vector<double> result = M.Gauss_Seidel(b, tolerance, x0);
        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        data.clear();
        x.clear();
        b.clear();
        i.clear();
        j.clear();
    }
}



TEST(Iteration_tests, Chebyshev_MPI_test) {
    std::size_t N = 100;
    std::size_t n = 10;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::ifstream fileL;
    std::ifstream fileL_max;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b;
    std::vector<double> x;
    std::vector<double> lambda_min_v;
    std::vector<double> lambda_max_v;
    fileA.open(py_path + "Iteration_tests/test_data.txt");
    filei.open(py_path + "Iteration_tests/test_i.txt");
    filej.open(py_path + "Iteration_tests/test_j.txt");
    fileb.open(py_path + "Iteration_tests/test_b.txt");
    filex.open(py_path + "Iteration_tests/test_x.txt");
    fileL.open(py_path + "Iteration_tests/test_l_min.txt");
    fileL_max.open(py_path + "Iteration_tests/test_L.txt");
    for (std::size_t it = 0; it < n; it++) {
        apply_vector<double>(fileA, data);
        apply_vector<double>(filex, x);
        apply_vector<double>(fileb, b);
        apply_vector<std::size_t>(filei, i);
        apply_vector<std::size_t>(filej, j);
        apply_vector<double>(fileL, lambda_min_v);
        apply_vector<double>(fileL_max, lambda_max_v);
        std::vector<DOK<double>> D;
        for (std::size_t z = 0; z < data.size(); z++) {
            D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        }
        CSR_matrix<double> M(D, N, N);
        double tolerance = pow(10, -10);
        std::vector<double> x0(N, 1);
        double lambda_r_tolerance0 = pow(10, -13);
        double lambda_max = M.estimate_lambda_max(x0, lambda_r_tolerance0);

        const std::size_t Number = 64;
        std::vector<double> result = M.Chebyshev_MPI(b, tolerance, x0, Number, lambda_max, lambda_min_v[0]);
        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }
        data.clear();
        x.clear();
        b.clear();
        i.clear();
        j.clear();
        lambda_max_v.clear();
        lambda_min_v.clear();
    }
}


TEST(Iteration_tests, SOR_test){
    std::size_t N = 100;
    std::size_t n = 10;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::ifstream fileL;
    std::ifstream fileW;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b;
    std::vector<double> x;
    std::vector<double> omega;
    fileA.open(py_path + "Chebyshev_SSOR_tests/test_data.txt");
    filei.open(py_path + "Chebyshev_SSOR_tests/test_i.txt");
    filej.open(py_path + "Chebyshev_SSOR_tests/test_j.txt");
    fileb.open(py_path + "Chebyshev_SSOR_tests/test_b.txt");
    filex.open(py_path + "Chebyshev_SSOR_tests/test_x.txt");
    fileL.open(py_path + "Chebyshev_SSOR_tests/test_L.txt");
    fileW.open(py_path + "Chebyshev_SSOR_tests/test_w.txt");

    for (std::size_t it = 0; it < n; it++) {
        apply_vector<double>(fileA, data);
        apply_vector<double>(filex, x);
        apply_vector<double>(fileb, b);
        apply_vector<std::size_t>(filei, i);
        apply_vector<std::size_t>(filej, j);
        apply_vector<double>(fileW, omega);
        std::vector<DOK<double>> D;
        for (std::size_t z = 0; z < data.size(); z++) {
            D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        }
        CSR_matrix<double> M(D, N, N);
        double tolerance = pow(10, -10);
        std::vector<double> x0(N);

        std::vector<double> result = M.SOR(b, tolerance, x0, omega[0]);
        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        data.clear();
        x.clear();
        b.clear();
        i.clear();
        j.clear();
        omega.clear();
    }
}



TEST(Iteration_tests, Chebyshev_Simmetrical_Gauss_Seidel_test){
    std::size_t N = 100;
    std::size_t n = 10;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::ifstream fileL;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b;
    std::vector<double> x;
    std::vector<double> L;
    fileA.open(py_path + "Iteration_Chebyshev_Gauss_Seidel_tests/test_data.txt");
    filei.open(py_path + "Iteration_Chebyshev_Gauss_Seidel_tests/test_i.txt");
    filej.open(py_path + "Iteration_Chebyshev_Gauss_Seidel_tests/test_j.txt");
    fileb.open(py_path + "Iteration_Chebyshev_Gauss_Seidel_tests/test_b.txt");
    filex.open(py_path + "Iteration_Chebyshev_Gauss_Seidel_tests/test_x.txt");
    fileL.open(py_path + "Iteration_Chebyshev_Gauss_Seidel_tests/test_L.txt");

    for (std::size_t it = 0; it < n; it++) {
        apply_vector<double>(fileA, data);
        apply_vector<double>(filex, x);
        apply_vector<double>(fileb, b);
        apply_vector<std::size_t>(filei, i);
        apply_vector<std::size_t>(filej, j);
        apply_vector<double>(fileL, L);
        std::vector<DOK<double>> D;
        for (std::size_t z = 0; z < data.size(); z++) {
            D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        }
        CSR_matrix<double> M(D, N, N);
        double tolerance = pow(10, -10);
        std::vector<double> x0(N);
        std::vector<double> result = M.Chebyshev_Simmetrical_Gauss_Seidel(b, tolerance, x0, L[0]);
        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        data.clear();
        x.clear();
        b.clear();
        i.clear();
        j.clear();
        L.clear();
    }
}


TEST(Iteration_tests, Chebyshev_SSOR_test){
    std::size_t N = 100;
    std::size_t n = 10;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::ifstream fileL;
    std::ifstream fileW;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b;
    std::vector<double> x;
    std::vector<double> spectr_radius;
    std::vector<double> omega;
    fileA.open(py_path + "Chebyshev_SSOR_tests/test_data.txt");
    filei.open(py_path + "Chebyshev_SSOR_tests/test_i.txt");
    filej.open(py_path + "Chebyshev_SSOR_tests/test_j.txt");
    fileb.open(py_path + "Chebyshev_SSOR_tests/test_b.txt");
    filex.open(py_path + "Chebyshev_SSOR_tests/test_x.txt");
    fileL.open(py_path + "Chebyshev_SSOR_tests/test_L.txt");
    fileW.open(py_path + "Chebyshev_SSOR_tests/test_w.txt");

    for (std::size_t it = 0; it < n; it++) {
        apply_vector<double>(fileA, data);
        apply_vector<double>(filex, x);
        apply_vector<double>(fileb, b);
        apply_vector<std::size_t>(filei, i);
        apply_vector<std::size_t>(filej, j);
        apply_vector<double>(fileL, spectr_radius);
        apply_vector<double>(fileW, omega);
        std::vector<DOK<double>> D;
        for (std::size_t z = 0; z < data.size(); z++) {
            D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        }
        CSR_matrix<double> M(D, N, N);
        double tolerance = pow(10, -10);
        std::vector<double> x0(N);

        std::vector<double> result = M.Chebyshev_SSOR(b, tolerance, x0, spectr_radius[0], omega[0]);
        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        data.clear();
        x.clear();
        b.clear();
        i.clear();
        j.clear();
        spectr_radius.clear();
        omega.clear();
    }
}

TEST(Iteration_tests, Steepest_Descent_test){
    std::size_t N = 100;
    std::size_t n = 10;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b;
    std::vector<double> x;
    fileA.open(py_path + "Chebyshev_SSOR_tests/test_data.txt");
    filei.open(py_path + "Chebyshev_SSOR_tests/test_i.txt");
    filej.open(py_path + "Chebyshev_SSOR_tests/test_j.txt");
    fileb.open(py_path + "Chebyshev_SSOR_tests/test_b.txt");
    filex.open(py_path + "Chebyshev_SSOR_tests/test_x.txt");
    for (std::size_t it = 0; it < n; it++) {
        apply_vector<double>(fileA, data);
        apply_vector<double>(filex, x);
        apply_vector<double>(fileb, b);
        apply_vector<std::size_t>(filei, i);
        apply_vector<std::size_t>(filej, j);

        std::vector<DOK<double>> D;
        for (std::size_t z = 0; z < data.size(); z++) {
            D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        }
        CSR_matrix<double> M(D, N, N);
        double tolerance = pow(10, -10);
        std::vector<double> x0(N);
        std::vector<double> result = M.Steepest_Descent(b, tolerance, x0);
        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        data.clear();
        x.clear();
        b.clear();
        i.clear();
        j.clear();
    }
}



TEST(Iteration_tests, Polak_Balls_test){
    std::size_t N = 100;
    std::size_t n = 10;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b;
    std::vector<double> x;
    fileA.open(py_path + "Chebyshev_SSOR_tests/test_data.txt");
    filei.open(py_path + "Chebyshev_SSOR_tests/test_i.txt");
    filej.open(py_path + "Chebyshev_SSOR_tests/test_j.txt");
    fileb.open(py_path + "Chebyshev_SSOR_tests/test_b.txt");
    filex.open(py_path + "Chebyshev_SSOR_tests/test_x.txt");
    for (std::size_t it = 0; it < n; it++) {
        apply_vector<double>(fileA, data);
        apply_vector<double>(filex, x);
        apply_vector<double>(fileb, b);
        apply_vector<std::size_t>(filei, i);
        apply_vector<std::size_t>(filej, j);

        std::vector<DOK<double>> D;
        for (std::size_t z = 0; z < data.size(); z++) {
            D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        }
        CSR_matrix<double> M(D, N, N);
        double tolerance = pow(10, -10);
        std::vector<double> x0(N, 1);
        std::vector<double> result = M.Polak_Balls(b, tolerance, x0);
        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        data.clear();
        x.clear();
        b.clear();
        i.clear();
        j.clear();
    }
}
TEST(Iteration_tests, Conjugate_Gradient_test){
    std::size_t N = 100;
    std::size_t n = 10;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b;
    std::vector<double> x;
    fileA.open(py_path + "Chebyshev_SSOR_tests/test_data.txt");
    filei.open(py_path + "Chebyshev_SSOR_tests/test_i.txt");
    filej.open(py_path + "Chebyshev_SSOR_tests/test_j.txt");
    fileb.open(py_path + "Chebyshev_SSOR_tests/test_b.txt");
    filex.open(py_path + "Chebyshev_SSOR_tests/test_x.txt");
    for (std::size_t it = 0; it < n; it++) {
        apply_vector<double>(fileA, data);
        apply_vector<double>(filex, x);
        apply_vector<double>(fileb, b);
        apply_vector<std::size_t>(filei, i);
        apply_vector<std::size_t>(filej, j);
        std::vector<DOK<double>> D;
        for (std::size_t z = 0; z < data.size(); z++) {
            D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        }
        CSR_matrix<double> M(D, N, N);
        double tolerance = pow(10, -10);
        std::vector<double> x0(N);

        std::vector<double> result = M.Conjugate_Gradient(b, tolerance, x0);
        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        data.clear();
        x.clear();
        b.clear();
        i.clear();
        j.clear();
    }
}


double average_time(auto&& data){
    double time = 0;
    for(auto&& it: data){
        time+=it.count();
    }
    return time/data.size();
}

TEST(Iteration_tests, time_test){
    std::size_t N = 100;
    std::size_t n = 10;
    std::ifstream fileA;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileb;
    std::ifstream filex;
    std::ifstream fileL_max;
    std::ifstream fileL_min;
    std::ifstream fileL_Cheb_GS;
    std::ifstream fileL_SYM_GS;
    std::ifstream fileL_SSOR;
    std::ifstream fileW;
    std::vector<double> data;
    std::vector<std::size_t> i;
    std::vector<std::size_t> j;
    std::vector<double> b;
    std::vector<double> x;
    std::vector<double> L_max;
    std::vector<double> L_min;
    std::vector<double> L_SSOR;
    std::vector<double> L_SYM_GS;
    std::vector<double> L_CHEB_GS;
    std::vector<double> omega;

    std::vector<std::chrono::duration<double>> Chebyshev_SSOR_time_v(n);
    std::vector<std::chrono::duration<double>> Conjugate_Gradient_time_v(n);
    std::vector<std::chrono::duration<double>> Polak_Balls_time_v(n);
    std::vector<std::chrono::duration<double>> Steepest_Descent_time_v(n);
    std::vector<std::chrono::duration<double>> Chebyshev_Simmetrical_Gauss_Seidel_time_v(n);
    std::vector<std::chrono::duration<double>> SOR_time_v(n);
    std::vector<std::chrono::duration<double>> MPI_time_v(n);
    std::vector<std::chrono::duration<double>> Jacobi_time_v(n);
    std::vector<std::chrono::duration<double>> Gauss_Seidel_time_v(n);
    std::vector<std::chrono::duration<double>> Chebyshev_MPI_time_v(n);

    fileA.open(py_path + "Universal_Iteraton_test/test_data.txt");
    filei.open(py_path + "Universal_Iteraton_test/test_i.txt");
    filej.open(py_path + "Universal_Iteraton_test/test_j.txt");
    fileb.open(py_path + "Universal_Iteraton_test/test_b.txt");
    filex.open(py_path + "Universal_Iteraton_test/test_x.txt");
    fileL_min.open(py_path + "Universal_Iteraton_test/test_L_min.txt");
    fileL_max.open(py_path + "Universal_Iteraton_test/test_L_max.txt");
    fileL_SSOR.open(py_path + "Universal_Iteraton_test/test_L_SSOR.txt");
    fileL_Cheb_GS.open(py_path + "Universal_Iteraton_test/test_L_CHEB_GS.txt");
    fileL_SYM_GS.open(py_path + "Universal_Iteraton_test/test_L_SYM_GS.txt");
    fileW.open(py_path + "Universal_Iteraton_test/test_w.txt");

    for (std::size_t it = 0; it < n; it++) {
        apply_vector<double>(fileA, data);
        apply_vector<double>(filex, x);
        apply_vector<double>(fileb, b);
        apply_vector<std::size_t>(filei, i);
        apply_vector<std::size_t>(filej, j);
        apply_vector<double>(fileL_SYM_GS, L_SYM_GS);
        apply_vector<double>(fileL_Cheb_GS, L_CHEB_GS);
        apply_vector<double>(fileL_max, L_max);
        apply_vector<double>(fileL_min, L_min);
        apply_vector<double>(fileL_SSOR, L_SSOR);
        apply_vector<double>(fileW, omega);
        std::vector<DOK<double>> D;
        for (std::size_t z = 0; z < data.size(); z++) {
            D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        }
        CSR_matrix<double> M(D, N, N);
        double tolerance = pow(10, -10);
        std::vector<double> x0(N, 1.0);

        auto start = std::chrono::steady_clock::now();
        std::vector<double> result = M.Chebyshev_SSOR(b, tolerance, x0, L_SSOR[0], omega[0]);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        Chebyshev_SSOR_time_v[it] = elapsed_seconds;

        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        double tau_MPI = 1/(L_min[0] + L_max[0]);
        start = std::chrono::steady_clock::now();
        result = M.MPI(b, tau_MPI, tolerance, x0);
        end = std::chrono::steady_clock::now();
        elapsed_seconds = end-start;
        MPI_time_v[it] = elapsed_seconds;

        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        start = std::chrono::steady_clock::now();
        result = M.Jacobi(b, tolerance,  x0);
        end = std::chrono::steady_clock::now();
        elapsed_seconds = end-start;
        Jacobi_time_v[it] = elapsed_seconds;

        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        start = std::chrono::steady_clock::now();
        result = M.Gauss_Seidel(b, tolerance, x0);
        end = std::chrono::steady_clock::now();
        elapsed_seconds = end-start;
        Gauss_Seidel_time_v[it] = elapsed_seconds;

        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        start = std::chrono::steady_clock::now();
        result = M.Chebyshev_Simmetrical_Gauss_Seidel(b, tolerance, x0, L_SYM_GS[0]);
        end = std::chrono::steady_clock::now();
        elapsed_seconds = end-start;
        Chebyshev_Simmetrical_Gauss_Seidel_time_v[it] = elapsed_seconds;//дописать

        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        start = std::chrono::steady_clock::now();
        result = M.Chebyshev_MPI(b, tolerance, x0, 64, L_max[0], L_min[0]); //переписать
        end = std::chrono::steady_clock::now();
        elapsed_seconds = end-start;
        Chebyshev_MPI_time_v[it] = elapsed_seconds;

        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        start = std::chrono::steady_clock::now();
        result = M.SOR(b, tolerance, x0, omega[0]); //переписать
        end = std::chrono::steady_clock::now();
        elapsed_seconds = end-start;
        SOR_time_v[it] = elapsed_seconds;

        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        start = std::chrono::steady_clock::now();
        result = M.Steepest_Descent(b, tolerance, x0);
        end = std::chrono::steady_clock::now();
        elapsed_seconds = end-start;
        Steepest_Descent_time_v[it] = elapsed_seconds;

        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        start = std::chrono::steady_clock::now();
        result = M.Polak_Balls(b, tolerance, x0);
        end = std::chrono::steady_clock::now();
        elapsed_seconds = end-start;
        Polak_Balls_time_v[it] = elapsed_seconds;

        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        start = std::chrono::steady_clock::now();
        result = M.Conjugate_Gradient(b, tolerance, x0);
        end = std::chrono::steady_clock::now();
        elapsed_seconds = end-start;
        Conjugate_Gradient_time_v[it] = elapsed_seconds;

        for (std::size_t p = 0; p < N; p++) {
            EXPECT_NEAR(result[p], x[p], tolerance);
        }

        data.clear();
        x.clear();
        b.clear();
        i.clear();
        j.clear();
        L_SYM_GS.clear();
        L_max.clear();
        L_min.clear();
        L_CHEB_GS.clear();
        L_SSOR.clear();
        omega.clear();
    }

    auto Chebyshev_SSOR_average = average_time(Chebyshev_SSOR_time_v);
    std::cout << "Chebyshev_SSOR: " << Chebyshev_SSOR_average<< "s\n";

    auto MPI_average = average_time(MPI_time_v);
    std::cout << "MPI: " << MPI_average<< "s\n";

    auto Jacobi_average = average_time(Jacobi_time_v);
    std::cout << "Jacobi: " << Jacobi_average<< "s\n";

    auto Gauss_Seidel_average = average_time(Gauss_Seidel_time_v);
    std::cout << "Gauss_Seidel: " << Gauss_Seidel_average<< "s\n";

    auto Chebyshev_MPI_average = average_time(Chebyshev_MPI_time_v);
    std::cout << "Chebyshev_MPI: " << Chebyshev_MPI_average<< "s\n";

    auto Chebyshev_Simmetrical_Gauss_Seidel_average = average_time(Chebyshev_Simmetrical_Gauss_Seidel_time_v);
    std::cout << "Chebyshev_Simmetrical_Gauss_Seidel: " << Chebyshev_Simmetrical_Gauss_Seidel_average<< "s\n";

    auto SOR_average = average_time(SOR_time_v);
    std::cout << "SOR: " << SOR_average<< "s\n";

    auto Steepest_Descent_average = average_time(Steepest_Descent_time_v);
    std::cout << "Steepest_Descent: " << Steepest_Descent_average<< "s\n";

    auto Polak_Balls_average = average_time(Polak_Balls_time_v);
    std::cout << "Polak_Balls: " << Polak_Balls_average<< "s\n";

    auto Conjugate_Gradient_average = average_time(Conjugate_Gradient_time_v);
    std::cout << "Conjugate_Gradient: " << Conjugate_Gradient_average<< "s\n";

}
