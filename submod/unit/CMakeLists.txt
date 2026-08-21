# file(GLOB_RECURSE unit_srcs src/*.cpp)
add_library(unit INTERFACE)
add_library(pemu::unit ALIAS unit)

find_package(mp-units CONFIG REQUIRED)
find_package(llnl-units CONFIG REQUIRED)

pemu_register_interface_headers(
    unit
    unit
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(
    unit
    INTERFACE
    llnl-units::units
    mp-units::mp-units
)

if(BUILD_TESTING)
    file(GLOB_RECURSE test_unit_srcs test/*.cpp)
    add_executable(test_unit ${test_unit_srcs})
    find_package(GTest CONFIG REQUIRED)
    target_link_libraries(test_unit PRIVATE GTest::gtest GTest::gtest_main unit)
    gtest_discover_tests(test_unit)
endif()
