file(GLOB_RECURSE mesh_srcs src/*.cpp)
add_library(mesh ${mesh_srcs})
add_library(pemu::mesh ALIAS mesh)
pemu_register_public_headers(
    mesh
    mesh
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
find_package(MOAB CONFIG REQUIRED)
target_link_libraries(mesh PUBLIC MOAB)

if(BUILD_TESTING)
    file(GLOB_RECURSE test_mesh_srcs test/*.cpp)
    add_executable(test_mesh ${test_mesh_srcs})
    find_package(GTest CONFIG REQUIRED)
    target_link_libraries(test_mesh PRIVATE GTest::gtest GTest::gtest_main mesh)
    target_compile_definitions(
        test_mesh
        PRIVATE
            PEMU_MESH_TEST_DATA_DIR="${PEMU_MESH_TEST_DATA_DIR}"  
    )
    gtest_discover_tests(test_mesh)
endif()