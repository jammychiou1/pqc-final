add_executable(mult_test mult_test.cpp)
#target_link_libraries(mult_test poc)
target_link_libraries(mult_test neon)
target_include_directories(mult_test PUBLIC ../include)

add_executable(mult_test_small mult_test_small.cpp)
#target_link_libraries(mult_test_small poc)
target_link_libraries(mult_test_small neon)
target_include_directories(mult_test_small PUBLIC ../include)

add_executable(ntt_10_profile ntt_10_profile.cpp)
target_link_libraries(ntt_10_profile neon)
target_include_directories(ntt_10_profile PUBLIC ../include)

add_executable(intt_10_x10_profile intt_10_x10_profile.cpp)
target_link_libraries(intt_10_x10_profile neon)
target_include_directories(intt_10_x10_profile PUBLIC ../include)

add_executable(ntt_9_profile ntt_9_profile.cpp)
target_link_libraries(ntt_9_profile neon)
target_include_directories(ntt_9_profile PUBLIC ../include)

add_executable(intt_9_x9_profile intt_9_x9_profile.cpp)
target_link_libraries(intt_9_x9_profile neon)
target_include_directories(intt_9_x9_profile PUBLIC ../include)
