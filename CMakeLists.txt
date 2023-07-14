cmake_minimum_required(VERSION 3.18)

project(Mult C CXX ASM)

option(OPTIMIZE "enable O3")

if(OPTIMIZE)
    add_compile_options(-O3)
endif()

option(DEBINFO "enable debug info")

if(OPTIMIZE)
    add_compile_options(-g)
endif()

if(IPO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT supported OUTPUT error)
    if( supported )
        message(STATUS "IPO / LTO enabled")
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(STATUS "IPO / LTO not supported: <${error}>")
    endif()
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(profiling)
