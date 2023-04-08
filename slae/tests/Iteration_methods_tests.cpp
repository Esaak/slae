#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include "../src/CSR_matrix.hpp"

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
    fileA.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_data.txt");
    filei.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_i.txt");
    filej.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_j.txt");
    fileb.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_b.txt");
    filex.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_x.txt");
    fileL.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_L.txt");

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
    fileA.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_data.txt");
    filei.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_i.txt");
    filej.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_j.txt");
    fileb.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_b.txt");
    filex.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_x.txt");

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
        result = M.Jacobi(b, data, tolerance);
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
    fileA.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_data.txt");
    filei.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_i.txt");
    filej.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_j.txt");
    fileb.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_b.txt");
    filex.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_x.txt");

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
    fileA.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_data.txt");
    filei.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_i.txt");
    filej.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_j.txt");
    fileb.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_b.txt");
    filex.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_x.txt");
    fileL.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_l_min.txt");
    fileL_max.open("/home/ilya/slae_lab/py_tests/Iteration_tests/test_L.txt");
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

        std::size_t Number = 64;
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
    fileA.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_data.txt");
    filei.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_i.txt");
    filej.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_j.txt");
    fileb.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_b.txt");
    filex.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_x.txt");
    fileL.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_L.txt");
    fileW.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_w.txt");

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
    fileA.open("/home/ilya/slae_lab/py_tests/Iteration_Chebyshev_Gauss_Seidel_tests/test_data.txt");
    filei.open("/home/ilya/slae_lab/py_tests/Iteration_Chebyshev_Gauss_Seidel_tests/test_i.txt");
    filej.open("/home/ilya/slae_lab/py_tests/Iteration_Chebyshev_Gauss_Seidel_tests/test_j.txt");
    fileb.open("/home/ilya/slae_lab/py_tests/Iteration_Chebyshev_Gauss_Seidel_tests/test_b.txt");
    filex.open("/home/ilya/slae_lab/py_tests/Iteration_Chebyshev_Gauss_Seidel_tests/test_x.txt");
    fileL.open("/home/ilya/slae_lab/py_tests/Iteration_Chebyshev_Gauss_Seidel_tests/test_L.txt");

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
    fileA.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_data.txt");
    filei.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_i.txt");
    filej.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_j.txt");
    fileb.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_b.txt");
    filex.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_x.txt");
    fileL.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_L.txt");
    fileW.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_w.txt");

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
    fileA.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_data.txt");
    filei.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_i.txt");
    filej.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_j.txt");
    fileb.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_b.txt");
    filex.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_x.txt");
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
    fileA.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_data.txt");
    filei.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_i.txt");
    filej.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_j.txt");
    fileb.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_b.txt");
    filex.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_x.txt");
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
    fileA.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_data.txt");
    filei.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_i.txt");
    filej.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_j.txt");
    fileb.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_b.txt");
    filex.open("/home/ilya/slae_lab/py_tests/Chebyshev_SSOR_tests/test_x.txt");
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

