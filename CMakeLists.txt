cmake_minimum_required(VERSION 3.30.0)

project(pemu VERSION 0.0.1 LANGUAGES C CXX)

if(BUILD_TESTING)
    include(CTest)
    enable_testing()
    include(GoogleTest)
endif()

add_subdirectory(submod)