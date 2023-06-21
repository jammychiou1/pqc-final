add_executable(mult_test mult_test.cpp)
target_link_libraries(mult_test poc)
target_include_directories(mult_test PUBLIC ../include)

add_executable(mult_test_small mult_test_small.cpp)
target_link_libraries(mult_test_small poc)
target_include_directories(mult_test_small PUBLIC ../include)
