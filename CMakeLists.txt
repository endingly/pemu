cmake_minimum_required(VERSION 3.30.0)

project(pemu VERSION 0.0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(BUILD_TESTING)
    include(CTest)
    enable_testing()
    include(GoogleTest)
endif()

add_subdirectory(submod)