# add_library(poc OBJECT poc/poc.cpp)
# target_include_directories(poc PUBLIC ../include)

add_library(neon OBJECT neon/neon.cpp neon/ntt_10.cpp neon/intt_10_x10.cpp neon/ntt_9.cpp neon/intt_9_x9.cpp neon/base_mul.cpp)
target_include_directories(neon PUBLIC ../include)
