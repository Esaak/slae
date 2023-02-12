# 'Google_test' is the subproject name

# 'lib' is the folder with Google Test sources
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_subdirectory(Google_tests/lib) # это строчка не нужна

add_executable(tridiagonal_test tridiagonal_test.cpp) # здесь название исполняемого файла и файл в котором лежит тест
add_test(NAME tridiagonal_test COMMAND tridiagonal_test.cpp)
target_link_libraries(tridiagonal_test slae gtest gtest_main)  # Здесь линковка, тот же исполняемый файл, название проекта и эти же либы

add_executable(CRS_matrix_test CRS_matrix_test.cpp)
add_test(NAME CRS_matrix_test COMMAND CRS_matrix_test.cpp)
target_link_libraries(CRS_matrix_test slae gtest gtest_main)


#set(i 0)
#file(GLOB files "*.cpp")
#foreach(filename ${files})
#    add_executable("test_${i}" ${filename})
#    add_test(NAME "test_${i}" COMMAND ${filename})
#    target_link_libraries("test_${i}" slae gtest gtest_main)
#    math(EXPR i "${i} + 1")
#endforeach()
