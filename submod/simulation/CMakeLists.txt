file(GLOB_RECURSE simulation_srcs src/*cpp)
add_library(simulation)
target_sources(simulation PRIVATE ${simulation_srcs})
add_library(pemu::simulation ALIAS simulation)
pemu_register_public_headers(
    simulation
    simulation
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(
    simulation
    PUBLIC
    pemu::equation
)

if(BUILD_TESTING)
    file(GLOB_RECURSE test_simulation_srcs test/*.cpp)
    add_executable(test_simulation ${test_simulation_srcs})
    find_package(GTest CONFIG REQUIRED)
    target_link_libraries(test_simulation PRIVATE GTest::gtest GTest::gtest_main simulation)
    target_compile_definitions(
        test_simulation
        PRIVATE
        PEMU_MESH_TEST_DATA_DIR="${PEMU_MESH_TEST_DATA_DIR}"
    )
    gtest_discover_tests(test_simulation)
endif()
