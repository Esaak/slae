#include <gtest/gtest.h>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include "../src/CSR_matrix.hpp"

using namespace DOK_space;
using namespace CSR_matrix_space;
const std::string py_path = "/home/ilya/SLAE/slae/py_tests/";
TEST(CSR_matrix_tests, DOK_sort){
    constexpr int N = 10;
    constexpr int minn = 1;
    constexpr int maxx = 1000;
    std::vector<int> i(N);
    std::vector<int> j(N);
    std::vector<int> data(N);
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<int>distrib(minn, maxx);
    std::vector<DOK<int>> D;
    for(auto& it: i){
        it = distrib(gen);
    }
    for(auto& it: j){
        it = distrib(gen);
    }
    for(auto& it: data){
        it = distrib(gen);
    }
    for(int z = 0; z<N; z++){
        D.emplace_back(DOK<int>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
    }
    std::sort(D.begin(), D.end());
    std::sort(i.begin(), i.end());
    std::vector<std::size_t> temp;
    std::transform(D.begin(), D.end(),std::back_inserter(temp), [](auto& a){return a.i;});

    ASSERT_EQ(temp.size(), i.size()) << "Vectors x and y are of unequal length";

    for (std::size_t p = 0; p < temp.size(); ++p) {
        EXPECT_EQ(i[p], temp[p]) << "Vectors x and y differ at index " << p;
    }

}

template<typename T>
void apply_vector(auto& file, std::vector<T>& data){
    std::string s, temp;
    getline(file, s);
    std::stringstream ss(s);
    while(ss>>temp){
        std::stringstream str;
        str<<temp;
        T f;
        str>>f;
        data.push_back(f);
    }
}

TEST(CSR_matrix_tests, boundary_condition){
    std::size_t N = 3;
    std::vector<std::size_t>indices{0, 1, 2};
    std::vector<std::size_t>indptr{0, 3, 3, 3};
    std::vector<std::size_t>i{0, 0, 0};
    std::vector<std::size_t>j{0, 1, 2};
    std::vector<double>data{1, 2, 3};
    std::vector<DOK<double>> D;
    for(std::size_t p = 0; p < N; p++){
        D.emplace_back(DOK<double>{static_cast<size_t>(i[p]), static_cast<size_t>(j[p]), data[p]});
    }
    CSR_matrix<double> M(D, N, N);
    for (std::size_t p = 0; p < N; p++) {
        EXPECT_DOUBLE_EQ(M(0, p), data[p]);
    }
    for (std::size_t p = 1; p < N; p++) {
        for(std::size_t q = 0; q < N; q++) {
            EXPECT_DOUBLE_EQ(M(p, q), 0);
        }
    }
    indptr.clear();
    indices.clear();
    i.clear();
    j.clear();
    data.clear();
    indices = {1};
    indptr = {0, 0, 1, 1};
    data = {1};
    i = {0};
    j = {0};
    std::vector<DOK<double>> D1;
    for(std::size_t p = 0; p < i.size(); p++){
        D1.emplace_back(DOK<double>{static_cast<size_t>(i[p]), static_cast<size_t>(j[p]), data[p]});
    }
    CSR_matrix<double> M1(D1, N, N);
    for(std::size_t p = 0; p < N; p++){
        for(std::size_t q = 0; q < N; q++){
            if(p == 0 && q == 0){
                EXPECT_DOUBLE_EQ(M1(p, q), 1);
            }
            else{
                EXPECT_DOUBLE_EQ(M1(p, q), 0);
            }
        }
    }

}
TEST(CSR_matrix_tests, index_operator_test){
    std::size_t N = 1000;
    std::size_t n = 10;
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
    fileA.open(py_path + "SRC_tests/test_data.txt");
    filed.open(py_path + "SRC_tests/test_indices.txt");
    filex.open(py_path + "SRC_tests/test_indptr.txt");
    filei.open(py_path + "SRC_tests/test_i.txt");
    filej.open(py_path + "SRC_tests/test_j.txt");
    for(std::size_t it = 0; it<n; it++) {
        apply_vector<double>(fileA, data);
        apply_vector<std::size_t>(filed, indices);
        apply_vector<std::size_t>(filex, indptr);
        apply_vector<std::size_t>(filei, i);
        apply_vector<std::size_t>(filej, j);
        std::vector<DOK<double>> D;
        for (std::size_t z = 0; z < data.size(); z++) {
            D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        }
        std::sort(D.begin(), D.end());
        CSR_matrix<double> M(D, N, N);

        for (std::size_t p = 0; p < data.size(); p++) {
            EXPECT_DOUBLE_EQ(M(D[p].i, D[p].j), D[p].value);
        }
        data.clear();
        indices.clear();
        indptr.clear();
        i.clear();
        j.clear();
        D.clear();
    }
    fileA.close();
    filed.close();
    filex.close();
    filei.close();
    filej.close();

}

TEST(CSR_matrix_tests, multiply_column){
    std::size_t n = 10;
    std::vector<double>data;
    std::vector<double>column;
    std::vector<double>X;
    std::vector<std::size_t>indices;
    std::vector<std::size_t>indptr;
    std::vector<std::size_t>i;
    std::vector<std::size_t>j;
    std::ifstream fileA;
    std::ifstream filed;
    std::ifstream filex;
    std::ifstream filei;
    std::ifstream filej;
    std::ifstream fileColumn;
    std::ifstream fileRes;
    fileA.open(py_path + "SRC_tests/test_data.txt");
    filed.open(py_path + "SRC_tests/test_indices.txt");
    filex.open(py_path + "SRC_tests/test_indptr.txt");
    filei.open(py_path + "SRC_tests/test_i.txt");
    filej.open(py_path + "SRC_tests/test_j.txt");
    fileColumn.open(py_path + "SRC_tests/test_D.txt");
    fileRes.open(py_path + "SRC_tests/test_X.txt");
    for(std::size_t it = 0; it< n; ++it) {
        apply_vector<double>(fileA, data);
        apply_vector<std::size_t>(filed, indices);
        apply_vector<std::size_t>(filex, indptr);
        apply_vector<std::size_t>(filei, i);
        apply_vector<std::size_t>(filej, j);
        apply_vector<double>(fileColumn, column);
        apply_vector<double>(fileRes, X);
        std::vector<DOK<double>> D;
        for (std::size_t z = 0; z < data.size(); z++) {
            D.emplace_back(DOK<double>{static_cast<size_t>(i[z]), static_cast<size_t>(j[z]), data[z]});
        }
        std::sort(D.begin(), D.end());
        CSR_matrix<double> M(D, X.size(), X.size());
        std::vector<double> Result = M * column;
        //EXPECT_EQ(Result, X);
        for (std::size_t p = 0; p < Result.size(); ++p) {
           EXPECT_NEAR(Result[p], X[p], std::abs(X[p])) << Result[p] << " " << X[p] << " " << p <<" "<<it<<"\n";
        }
        i.clear();
        j.clear();
        indptr.clear();
        indices.clear();
        X.clear();
        column.clear();
        data.clear();
    }
    fileA.close();
    filed.close();
    filex.close();
    filei.close();
    filej.close();
    fileColumn.close();
    fileRes.close();
}

TEST(CSR_matrix_tests, transpose_test){
    std::size_t N = 3;
    std::vector<std::size_t>indices{0, 1, 2};
    std::vector<std::size_t>indptr{0, 3, 3, 3};
    std::vector<std::size_t>i{0, 0, 0};
    std::vector<std::size_t>j{0, 1, 2};
    std::vector<double>data{1, 2, 3};
    std::vector<DOK<double>> D;
    for(std::size_t p = 0; p < N; p++){
        D.emplace_back(DOK<double>{static_cast<size_t>(i[p]), static_cast<size_t>(j[p]), data[p]});
    }
    CSR_matrix<double> M(D, N, N);
    CSR_matrix<double> B = M.transpose();
    indptr.clear();
    indices.clear();
    i.clear();
    j.clear();
    data.clear();
    N = 4;

    indices = {0,2,0,3,0,1,2,3};
    indptr = {0,2,2,4,8};
    data = {1,4,5,8,2,3,5,6};
    std::vector<double> dataT = {1,5,2,3,4,5,8,6};

    i = {0, 0, 2, 2, 3, 3, 3, 3};
    std::vector<double> iT = {0, 0, 0, 1, 2, 2, 3, 3};
    j = {0, 2, 0, 3, 0, 1, 2, 3};
    std::vector<double> jT = {0, 2, 3, 3, 0, 3, 2, 3};
    std::vector<DOK<double>> D1;
    for(std::size_t p = 0; p < i.size(); p++){
        D1.emplace_back(DOK<double>{static_cast<size_t>(i[p]), static_cast<size_t>(j[p]), data[p]});
    }
    CSR_matrix<double> M1(D1, N, N);
    CSR_matrix<double> B1 = M1.transpose();
    for(std::size_t it = 0; it < iT.size(); it++){
        EXPECT_DOUBLE_EQ(B1(iT[it], jT[it]), dataT[it]);
    }

}