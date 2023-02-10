cmake_minimum_required(VERSION 3.24)
project(slae)


set(CMAKE_CXX_STANDARD 23)
add_library(slae INTERFACE tridiagonal_matrix.hpp slae_solver.hpp )

