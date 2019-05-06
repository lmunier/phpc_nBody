cmake_minimum_required(VERSION 3.13)
project(Project)

set(CMAKE_CXX_STANDARD 17)

add_executable(Project main.cpp Vector.hpp Particle.cpp Particle.hpp constants.hpp)