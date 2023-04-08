cmake_minimum_required(VERSION 3.24)
project(slae)


set(CMAKE_CXX_STANDARD 23)
add_compile_options(-Wall -Werror -O3)
add_subdirectory(tests)
add_subdirectory(src)
