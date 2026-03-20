cmake_minimum_required(VERSION 4.0)
project(assgn)

set(CMAKE_CXX_STANDARD 20)

add_executable(assgn main.cpp
        Bhtree.cpp Octant.cpp Bhtree.h Collisions.h Constants.h Octant.h Types.h)
