# file(GLOB_RECURSE bc_srcs src/*cpp)
add_library(boundary INTERFACE)
add_library(pemu::boundary ALIAS boundary)
pemu_register_interface_headers(
    boundary
    boundary
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(
    boundary INTERFACE pemu::mesh
)

if(BUILD_TESTING)
    file(GLOB_RECURSE test_boundary_srcs test/*.cpp)
    add_executable(test_boundary ${test_boundary_srcs})
    find_package(GTest CONFIG REQUIRED)
    target_link_libraries(test_boundary PRIVATE GTest::gtest GTest::gtest_main boundary)
    target_compile_definitions(
        test_boundary
        PRIVATE
        PEMU_MESH_TEST_DATA_DIR="${PEMU_MESH_TEST_DATA_DIR}"
    )
    gtest_discover_tests(test_boundary)
endif()
