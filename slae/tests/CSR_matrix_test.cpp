#include <gtest/gtest.h>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include "../src/CSR_matrix.hpp"

using namespace DOK_space;
using namespace CSR_matrix_space;

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

TEST(CSR_matrix_tests, index_operator_test){
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
    filej.open("/home/ilya/slae_lab/py_tests/SRC_tests/test_j.txt");
    for(std::size_t it; it<N; it++) {
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
    std::size_t N = 10;
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
    fileA.open("/home/ilya/slae_lab/py_tests/SRC_tests/test_data.txt");
    filed.open("/home/ilya/slae_lab/py_tests/SRC_tests/test_indices.txt");
    filex.open("/home/ilya/slae_lab/py_tests/SRC_tests/test_indptr.txt");
    filei.open("/home/ilya/slae_lab/py_tests/SRC_tests/test_i.txt");
    filej.open("/home/ilya/slae_lab/py_tests/SRC_tests/test_j.txt");
    fileColumn.open("/home/ilya/slae_lab/py_tests/SRC_tests/test_D.txt");
    fileRes.open("/home/ilya/slae_lab/py_tests/SRC_tests/test_X.txt");
    for(std::size_t it = 0; it< N; ++it) {
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
        CSR_matrix<double> M(D, N, N);
        std::vector<double> Result = M * column;
        //EXPECT_EQ(Result, X);
        for (std::size_t p = 0; p < Result.size(); ++p) {
           EXPECT_DOUBLE_EQ(Result[p], X[p]) << Result[p] << " " << X[p] << " " << p <<" "<<it<<"\n";
        }
        D.clear();
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
