# 'Google_test' is the subproject name

# 'lib' is the folder with Google Test sources
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_subdirectory(Google_tests/lib) # это строчка не нужна


file(GLOB files "*.cpp" "./TEST_04_03/*.cpp" "./TEST_08_04/*.cpp")
#file(GLOB files1 "./TEST_04_03/*.cpp")
foreach(filename ${files})
    get_filename_component (TName ${filename} NAME_WE)
    add_executable("${TName}" ${filename})
    add_test(NAME "${TName}" COMMAND ${filename})
    target_link_libraries(${TName} slae gtest gtest_main)
endforeach()

