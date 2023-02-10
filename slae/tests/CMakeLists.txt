# 'Google_test' is the subproject name

# 'lib' is the folder with Google Test sources
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_subdirectory(Google_tests/lib)

add_executable(tridiagonal_test tridiagonal_test.cpp)
add_test(NAME tridiagonal_test COMMAND tridiagonal_test.cpp)
target_link_libraries(tridiagonal_test slae gtest gtest_main)

add_executable(test test_tesT.cpp)
add_test(NAME test COMMAND test_tesT.cpp)
target_link_libraries(test slae gtest gtest_main)

# foreach
# regax in cmake
# флаги cmake (-D{название флага})
