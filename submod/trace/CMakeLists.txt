# file(GLOB_RECURSE trace_srcs src/*.cpp)
add_library(trace INTERFACE)
add_library(pemu::trace ALIAS trace)

pemu_register_interface_headers(
    trace
    trace
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

if(BUILD_TESTING)
    file(GLOB_RECURSE test_trace_srcs test/*.cpp)
    add_executable(test_trace ${test_trace_srcs})
    find_package(GTest CONFIG REQUIRED)
    target_link_libraries(test_trace PRIVATE GTest::gtest GTest::gtest_main trace)
    gtest_discover_tests(test_trace)
endif()
