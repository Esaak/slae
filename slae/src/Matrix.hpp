
#ifndef SLAE_MATRIX_HPP


#define SLAE_MATRIX_HPP


#include <vector>
#include <functional>
#include <numeric>
#include <utility>
#include <cmath>
#include <compare>
#include <cstring>
#include <algorithm>


template<typename T>
using IsArithmetical = std::enable_if_t<std::is_arithmetic_v<T>, bool>;

namespace Mrx {
    template<typename T, IsArithmetical<T> = true>
    class Matrix {
    private:
        std::size_t column{};
        std::size_t row{};
        std::vector<T> data;
    public:
        Matrix() = default;
        explicit Matrix(const std::vector<std::vector<T>> &other):column(other[0].size()), row(other.size()){
            for (std::size_t i = 0; i < other.size(); ++i) {
                for (std::size_t j = 0; j < other[0].size(); ++j) {
                    data.push_back(other[i][j]);
                }
            }
        }


        T &operator()(std::size_t i, std::size_t j) {
            return data[i * column + j];
        }

        T operator()(std::size_t i, std::size_t j) const {
            return data[i * column + j];
        }

        Matrix &operator*=(const Matrix<T> &other) {
            if (column != other.row) throw std::invalid_argument("Invalid");
            std::vector<T> new_data(row * other.column);
            for (std::size_t i = 0; i < row; i++) {
                for (std::size_t j = 0; j < other.column; j++) {
                    for (std::size_t k = 0; k < column; k++) {
                        new_data[j + i * other.column] += data[k + i * column] * other.data[j + k * other.column];
                    }
                }
            }
            column = other.column;
            data = new_data;
            return *this;
        }

        const std::size_t get_column_size() const {
            return column;
        }

        const std::size_t get_row_size() const {
            return row;
        }


        Matrix<T> transponse() const{
            Matrix<T> new_matrix;
            new_matrix.data.resize(row * column);
            new_matrix.column = row;
            new_matrix.row = column;
            if (column == row) {
                for (std::size_t i = 0; i < row; i++) {
                    for (std::size_t j = 1 + i; j < column; j++) {
                        new_matrix(j, i) = data[i * column + j];
                        new_matrix(i, j) = data[j * column + i];
                    }
                    new_matrix(i, i) = data[i * column + i];
                }
                return new_matrix;
            }
            for (std::size_t t = 0; t < column; t++) {
                for (std::size_t q = 0; q < row; q++) {
                    new_matrix(t, q) = data[q * column + t];
                }
            }
            return new_matrix;
        }

        Matrix &operator*(const Matrix<T> &other) const{
            if (column != other.row) throw std::invalid_argument("Invalid");
            Matrix new_matrix = *this;
            return new_matrix *= other;
        }

        std::vector<T> dot(const std::vector<T> &vec) const{
            if (column != vec.size()) throw std::invalid_argument("Invalid");
            std::vector<T> res_vector(row);
            for (std::size_t i = 0; i < vec.size(); ++i) {
                res_vector[i] = std::inner_product(data.begin() + column * i, data.begin() + column * (i + 1),
                                                   vec.begin(), static_cast<T>(0));
            }
            return res_vector;
        }

        Matrix &operator+=(const Matrix<T> &other) {
            data += other.data;
            return *this;
        }

        Matrix operator+(const Matrix<T> &other) const{
            Matrix res_matrix = *this;
            return res_matrix += other;
        }

        std::vector<T> operator[](std::size_t i) const {
            std::vector<T> one_row;
            std::copy_n(data.begin() + column * i, column, std::back_inserter(one_row));
            return one_row;
        }

        std::vector<T> get_row(std::size_t i, int begin, int end) const {
            std::vector<T> x(end - begin);
            std::copy(data.begin() + i * column + begin, data.begin() + i * column + end, x.begin());
            return x;
        }

        std::vector<T> get_column(std::size_t j, int begin, int end) const {
            std::vector<T> x(end - begin);
            for (int i = begin; i < end; i++) {
                x[i - begin] = (*this)(i, j);
            }
            return x;
        }

        static Matrix eye(std::size_t i){
            Matrix<T> new_matrix;
            new_matrix.data.resize(i*i);
            new_matrix.row = i;
            new_matrix.column = i;
            for (std::size_t p = 0; p < i; p++) {
                new_matrix(p,p) = 1;
            }
            return new_matrix;
        }

    };


}
#endif //SLAE_MATRIX_HPP
