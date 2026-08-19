# file(GLOB_RECURSE physics_srcs src/*cpp)
add_library(physics INTERFACE)
# target_sources(physics PRIVATE ${physics_srcs})
add_library(pemu::physics ALIAS physics)
pemu_register_interface_headers(
    physics
    physics
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(
    physics
    INTERFACE
    pemu::field
)

if(BUILD_TESTING)
    file(GLOB_RECURSE test_physics_srcs test/*.cpp)
    add_executable(test_physics ${test_physics_srcs})
    find_package(GTest CONFIG REQUIRED)
    target_link_libraries(test_physics PRIVATE GTest::gtest GTest::gtest_main physics)
    target_compile_definitions(
        test_physics
        PRIVATE
        PEMU_MESH_TEST_DATA_DIR="${PEMU_MESH_TEST_DATA_DIR}"
    )
    gtest_discover_tests(test_physics)
endif()
