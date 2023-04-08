# CMake generated Testfile for 
# Source directory: /home/ilya/slae_lab/slae/tests
# Build directory: /home/ilya/slae_lab/slae/cmake-build-debug/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[test_0]=] "/home/ilya/slae_lab/slae/tests/CRS_matrix_test.cpp")
set_tests_properties([=[test_0]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/ilya/slae_lab/slae/tests/CMakeLists.txt;19;add_test;/home/ilya/slae_lab/slae/tests/CMakeLists.txt;0;")
add_test([=[test_1]=] "/home/ilya/slae_lab/slae/tests/tridiagonal_test.cpp")
set_tests_properties([=[test_1]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/ilya/slae_lab/slae/tests/CMakeLists.txt;19;add_test;/home/ilya/slae_lab/slae/tests/CMakeLists.txt;0;")
subdirs("Google_tests/lib")
