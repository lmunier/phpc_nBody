cmake_minimum_required(VERSION 3.14)
project(phpc_nBody)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

add_executable(phpc_nBody
        constants.hpp
        main.cpp
        Tree.hpp
        Tree.cpp
        Vector.hpp)
