
#ifndef SLAE_MATRIX_HPP


#define SLAE_MATRIX_HPP

#include <concepts>
#include <ranges>
#include <vector>
#include <functional>
#include <numeric>
#include <utility>
#include <cmath>
template<typename T>
concept arithmetical = std::is_floating_point<T>::value || std::is_integral<T>::value;


namespace Mrx {
    template<arithmetical T>
    class Matrix {
    private:
        std::size_t column;
        std::size_t row;
        std::vector<T> data;

    public:
        void change_matrix(const std::vector<std::vector<T>> &other) {
            data.clear();
            for (std::size_t i = 0; i < other.size(); ++i) {
                for (std::size_t j = 0; j < other[0].size(); ++j) {
                    data.push_back(other[i][j]);
                }
            }
            row = other.size();
            column = other[0].size();
        }

        T& operator()(std::size_t i, std::size_t j){
            return data[i * column + j];
        }
        T operator()(std::size_t i, std::size_t j) const{
            return data[i*column  + j];
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
        const std::size_t get_column_size() const{
            return column;
        }
        const std::size_t get_row_size() const{
            return row;
        }

        //friend Matrix<T> operator * (const Matrix<T>& A, const Matrix<T>& B);
        void transponse(){
            std::vector<T>new_data(row*column);
            if(column == row){
            for(std::size_t i = 0; i < row; i++){
                for(std::size_t j = 1+i; j < column; j++){
                    T tmp = data[i * row + j];
                    new_data[i*row+j] = data[j*row+i];
                    new_data[j*row+i] = tmp;
                }
            }
            for(std::size_t i=0; i<column; i++){
                new_data[i*column + i]=data[i*column + i];
            }
            data = new_data;
            return;
        }
        for(std::size_t t=0; t<column; t++){
            for(std::size_t q=0;q<row; q++){
                new_data[t*row + q] =data[q*column + t];
            }
        }
        unsigned temp= row;
        row = column;
        column = temp;
        data = new_data;
        }

        Matrix &operator * (const Matrix<T>& other){
            if (column != other.row) throw std::invalid_argument("Invalid");
            Matrix new_matrix = *this;
            return new_matrix*=other;
        }
        std::vector<T> dot(const std::vector<T>& vec){
            if(column != vec.size()) throw std::invalid_argument("Invalid");
            std::vector<T> res_vector(row);
            for(std::size_t i = 0; i < vec.size(); ++i){
                res_vector[i] = std::inner_product(data.begin() + column * i, data.begin() + column * (i + 1), vec.begin(), static_cast<T>(0));
            }
            return res_vector;
        }

        Matrix& operator += (const Matrix<T>& other){
            if (column != other.column || row != other.row) throw std::invalid_argument("Invalid");
            data+=other.data;
            return *this;
        }
        Matrix operator + (const Matrix<T>& other){
            Matrix res_matrix = *this;
            return res_matrix+= other;
        }

        std::vector<T> operator [](std::size_t i) const{
            std::vector<T> one_row;
            std::ranges::copy_n(data.begin() + column*i, column, std::back_inserter(one_row));
            return one_row;
        }

        std::vector<T> get_row(std:: size_t i, int begin, int end) const{
            std::vector<T> x(end - begin);
            std::ranges::copy(data.begin() + i*column + begin, data.begin() + i*column + end, x.begin());
            return x;
        }

        std::vector<T> get_column(std:: size_t j, int begin, int end) const{
            std::vector<T> x(end - begin);
            for(int i = begin; i < end; i++){
                x[i - begin]=(*this)(i,j);
            }
            return x;
        }
        void ones(std::size_t j, std::size_t i){
            data.resize(i*j);
            row = i;
            column = j;
            for(std::size_t p = 0; p < i; p++){
                data[p + p*j] = 1;
            }
        }

    };
    /*
    template<aritmetical T>
    Matrix<T> operator * (const Matrix<T>& A, const Matrix<T>& B){
        Matrix<T> new_matrix = A;
        return new_matrix*=B;
    }*/

}
#endif //SLAE_MATRIX_HPP
