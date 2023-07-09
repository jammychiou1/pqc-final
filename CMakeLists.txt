cmake_minimum_required(VERSION 3.18)

project(Mult C CXX ASM)

# option(OPTIMIZE "enable O3")
#
# if(OPTIMIZE)
#     add_compile_options(-O3)
# endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(profiling)
