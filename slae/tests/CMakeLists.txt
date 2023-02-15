# 'Google_test' is the subproject name

# 'lib' is the folder with Google Test sources
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_subdirectory(Google_tests/lib) # это строчка не нужна

#add_executable(tridiagonal_test tridiagonal_test.cpp) # здесь название исполняемого файла и файл в котором лежит тест
#add_test(NAME tridiagonal_test COMMAND tridiagonal_test.cpp)
#target_link_libraries(tridiagonal_test slae gtest gtest_main)  # Здесь линковка, тот же исполняемый файл, название проекта и эти же либы
#
#add_executable(CRS_matrix_test CRS_matrix_test.cpp)
#add_test(NAME CRS_matrix_test COMMAND CRS_matrix_test.cpp)
#target_link_libraries(CRS_matrix_test slae gtest gtest_main)


file(GLOB files "*.cpp")
foreach(filename ${files})
    get_filename_component (TName ${filename} NAME_WE)
    add_executable("${TName}" ${filename})
    add_test(NAME "${TName}" COMMAND ${filename})
    target_link_libraries(${TName} slae gtest gtest_main)
endforeach()

#create_test_sourcelist (Tests
#        CommonCxxTests.cxx
#        CRS_matrix_test.cpp
#        tridiagonal_test.cpp
#        )
#
## add the executable
#add_executable (CommonCxxTests ${Tests})
#
## remove the test driver source file
#set (TestsToRun ${Tests})
#remove (TestsToRun CommonCxxTests.cxx)
#
## Add all the ADD_TEST for each test
#foreach (test ${TestsToRun})
#    get_filename_component (TName ${test} NAME_WE)
#    add_test (NAME ${TName} COMMAND CommonCxxTests ${TName})
#endforeach ()