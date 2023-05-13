#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include "../src/CSR_matrix.hpp"

using namespace DOK_space;
using namespace CSR_matrix_space;

TEST(Iteration_simple_tests, lambda_max_simple_test) {
    std::size_t N = 9;
    std::vector<std::size_t> i{0, 0, 0, 1, 1, 1, 2, 2, 2};
    std::vector<std::size_t> j{0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> data{200, 1, 2, 3, 200, 4, 5, 6, 200};
    std::vector<DOK<double>> D;
    for (std::size_t p = 0; p < N; p++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[p]), static_cast<size_t>(j[p]), data[p]});
    }
    std::vector<double> x0(3, 1);
    CSR_matrix<double> M(D, 3, 3);
    double lambda_r_tolerance0 = pow(10, -6);
    double lambda_max = M.estimate_lambda_max(x0, lambda_r_tolerance0);
    std::cout << lambda_max;
}

TEST(Iteration_simple_tests, Chebyshev_MPI_simple_test) {
    std::vector<DOK<double>> D;
    std::vector<int> i{0, 0, 0, 1, 1, 1, 2, 2, 2};
    std::vector<int> j{0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> data{12, 17, 3, 17, 15825, 28, 3, 28, 238};
    std::vector<double> x0{1, 1, 1};
    std::vector<double> b{1, 2, 3};
    for (std::size_t z = 0; z < i.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    CSR_matrix<double> ANSW(D, 3, 3);
    double lambda_min = 11.9427;
    double lambda_r_tolerance0 = pow(10, -13);
    double lambda_max = ANSW.estimate_lambda_max(x0, lambda_r_tolerance0);

    std::vector<double> a = ANSW.Chebyshev_MPI(b, lambda_r_tolerance0, x0, 128, lambda_max, lambda_min);
    std::cout << "\n";
    std::copy(a.begin(), a.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";

}

TEST(Iteration_simple_tests, SOR_simple_test) {
    std::vector<DOK<double>> D;
    std::vector<int> i{0, 0, 0, 1, 1, 1, 2, 2, 2};
    std::vector<int> j{0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> data{12, 17, 3, 17, 15825, 28, 3, 28, 238};
    std::vector<double> x0{1, 1, 1};
    std::vector<double> b{1, 2, 3};
    std::vector<double> diag_elements(3);
    for (std::size_t z = 0; z < i.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        if (i[z] == j[z]) {
            diag_elements[i[z]] = data[z];
        }
    }
    CSR_matrix<double> ANSW(D, 3, 3);
    double lambda_r_tolerance0 = pow(10, -13);
    for (auto &p: D) {
        if (p.i == p.j) p.value = 0;
        else {
            p.value /= -diag_elements[p.i];
        }
    }
    CSR_matrix<double> C(D, 3, 3);
    double lambda_max = C.estimate_lambda_max(x0, lambda_r_tolerance0);
    double omega = 1 + std::pow(lambda_max / (1 + std::sqrt(1 - std::pow(lambda_max, 2))), 2);

    std::vector<double> a = ANSW.SOR(b, lambda_r_tolerance0, x0, omega);
    std::cout << "\n";
    std::copy(a.begin(), a.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
}

TEST(Iteration_simple_tests, Chebyshev_Simmetrical_Gauss_Seidel_simple_test) {
    std::vector<DOK<double>> D;
    std::vector<int> i{0, 0, 0, 1, 1, 1, 2, 2, 2};
    std::vector<int> j{0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> data{6, 2, 1, 3, 7, 4, 1, 5, 8};
    std::vector<double> x0{1, 1, 1};
    std::vector<double> b{1, 2, 3};
    for (std::size_t z = 0; z < i.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    CSR_matrix<double> ANSW(D, 3, 3);
    double tolerance = std::pow(10, -13);
    double spectr_radius = 0.2336;
    std::vector<double> a = ANSW.Chebyshev_Simmetrical_Gauss_Seidel(b, tolerance, x0, spectr_radius);
    std::cout << "\n";
    std::copy(a.begin(), a.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
}

TEST(Iteration_simple_tests, Steeptest_Descent_simple_test) {
    std::vector<DOK<double>> D;
    std::vector<int> i{0, 0, 0, 1, 1, 1, 2, 2, 2};
    std::vector<int> j{0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> data{12, 17, 3, 17, 15825, 28, 3, 28, 238};
    std::vector<double> x0{1, 1, 1};
    std::vector<double> b{1, 2, 3};
    std::vector<double> diag_elements(3);
    for (std::size_t z = 0; z < i.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        if (i[z] == j[z]) {
            diag_elements[i[z]] = data[z];
        }
    }
    CSR_matrix<double> ANSW(D, 3, 3);
    double lambda_r_tolerance0 = pow(10, -13);
    std::vector<double> a = ANSW.Steepest_Descent(b, lambda_r_tolerance0, x0);
    std::cout << "\n";
    std::copy(a.begin(), a.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
}

TEST(Iteration_simple_tests, Polak_Balls_simple_test) {
    std::vector<DOK<double>> D;
    std::vector<int> i{0, 0, 0, 1, 1, 1, 2, 2, 2};
    std::vector<int> j{0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> data{12, 17, 3, 17, 15825, 28, 3, 28, 238};
    std::vector<double> x0{0.01, 0.01, 0.01};
    std::vector<double> b{1, 2, 3};
    std::vector<double> diag_elements(3);
    for (std::size_t z = 0; z < i.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        if (i[z] == j[z]) {
            diag_elements[i[z]] = data[z];
        }
    }
    CSR_matrix<double> ANSW(D, 3, 3);
    double lambda_r_tolerance0 = pow(10, -13);
    std::vector<double> a = ANSW.Polak_Balls(b, lambda_r_tolerance0, x0);
    std::cout << "\n";
    std::copy(a.begin(), a.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
}

TEST(Iteration_simple_tests, Conjugate_Gradient_simple_test) {
    std::vector<DOK<double>> D;
    std::vector<int> i{0, 0, 0, 1, 1, 1, 2, 2, 2};
    std::vector<int> j{0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> data{12, 17, 3, 17, 15825, 28, 3, 28, 238};
    std::vector<double> x0{1, 1, 1};
    std::vector<double> b{1, 2, 3};
    std::vector<double> diag_elements(3);
    for (std::size_t z = 0; z < i.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        if (i[z] == j[z]) {
            diag_elements[i[z]] = data[z];
        }
    }
    CSR_matrix<double> ANSW(D, 3, 3);
    double lambda_r_tolerance0 = pow(10, -13);
    std::vector<double> a = ANSW.Conjugate_Gradient(b, lambda_r_tolerance0, x0);
    std::cout << "\n";
    std::copy(a.begin(), a.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
}


TEST(Iteration_simple_tests, BCG_simple_test) {
    std::vector<DOK<double>> D;
    std::vector<int> i{0, 0, 0, 1, 1, 1, 2, 2, 2};
    std::vector<int> j{0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> data{12, 17, 3, 17, 15825, 28, 3, 28, 238};
    std::vector<double> x0{0, 0, 0};
    std::vector<double> b{1, 2, 3};
    std::vector<double> diag_elements(3);
    for (std::size_t z = 0; z < i.size(); z++) {
        D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        if (i[z] == j[z]) {
            diag_elements[i[z]] = data[z];
        }
    }
    CSR_matrix<double> ANSW(D, 3, 3);
    double lambda_r_tolerance0 = pow(10, -13);
    std::vector<double> a = ANSW.BCG(b, lambda_r_tolerance0, x0);
    std::cout << "\n";
    std::copy(a.begin(), a.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
}
