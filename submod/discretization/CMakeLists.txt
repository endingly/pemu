file(GLOB_RECURSE discretization_srcs src/*cpp)
add_library(discretization)
target_sources(discretization PRIVATE ${discretization_srcs})
add_library(pemu::discretization ALIAS discretization)
pemu_register_public_headers(
    discretization
    discretization
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(
    discretization
    PUBLIC
    pemu::boundary
    pemu::field
    pemu::linalg
    pemu::mesh
    pemu::physics
)

if(BUILD_TESTING)
    file(GLOB_RECURSE test_discretization_srcs test/*.cpp)
    add_executable(test_discretization ${test_discretization_srcs})
    find_package(GTest CONFIG REQUIRED)
    target_link_libraries(test_discretization PRIVATE GTest::gtest GTest::gtest_main discretization)
    target_compile_definitions(
        test_discretization
        PRIVATE
        PEMU_MESH_TEST_DATA_DIR="${PEMU_MESH_TEST_DATA_DIR}"
    )
    gtest_discover_tests(test_discretization)
endif()
