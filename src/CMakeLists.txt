# add_library(poc OBJECT poc/poc.cpp)
# target_include_directories(poc PUBLIC ../include)

add_library(neon OBJECT neon/neon.cpp neon/ntt_10.cpp neon/intt_10_x10.cpp neon/ntt_9.cpp neon/intt_9_x9.cpp neon/base_mul.cpp neon/mult_low.cpp)
target_include_directories(neon PUBLIC ../include)

add_library(rader_outer OBJECT rader_outer/polymul.c rader_outer/ntt_17.c rader_outer/intt_17.c rader_outer/rader.S rader_outer/ntt_3_2.c rader_outer/intt_3_2.c rader_outer/radix3.S rader_outer/base_mul.c rader_outer/mul.S)
target_include_directories(rader_outer PUBLIC ../include)
