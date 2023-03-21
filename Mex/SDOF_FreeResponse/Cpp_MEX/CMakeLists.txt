cmake_minimum_required(VERSION 3.0.0)
project(CPP_MEX VERSION 0.1.0 LANGUAGES CXX)

file(GLOB_RECURSE SOURCE_FILES FOLLOW_SYMLINKS
    ./*.cpp
)

include_directories("$(MATLABROOT)/extern/include")

add_library(CPP_MEX SHARED ${SOURCE_FILES})
