add_library(poc OBJECT poc/poc.cpp)
target_include_directories(poc PUBLIC ../include)

add_library(neon OBJECT neon/neon.cpp)
target_include_directories(neon PUBLIC ../include)
