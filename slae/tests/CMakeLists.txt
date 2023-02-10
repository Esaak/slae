# 'Google_test' is the subproject name

# 'lib' is the folder with Google Test sources
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_subdirectory(Google_tests/lib) # это строчка не нужна

add_executable(tridiagonal_test tridiagonal_test.cpp) # здесь название исполняемого файла и файл в котором лежит тест
add_test(NAME tridiagonal_test COMMAND tridiagonal_test.cpp)
target_link_libraries(tridiagonal_test slae gtest gtest_main)  # Здесь линковка, тот же исполняемый файл, название проекта и эти же либы



# foreach
# regax in cmake
# флаги cmake (-D{название флага})
