add_executable(mult_test mult_test.cpp)
#target_link_libraries(mult_test poc)
target_link_libraries(mult_test neon)
target_include_directories(mult_test PUBLIC ../include)

add_executable(mult_test_small mult_test_small.cpp)
#target_link_libraries(mult_test_small poc)
target_link_libraries(mult_test_small neon)
target_include_directories(mult_test_small PUBLIC ../include)
