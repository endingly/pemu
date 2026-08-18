# file(GLOB_RECURSE field_srcs src/*cpp)
add_library(field INTERFACE)
add_library(pemu::field ALIAS field)
pemu_register_interface_headers(
    field
    field
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(
    field INTERFACE pemu::mesh
)

if(BUILD_TESTING)
    file(GLOB_RECURSE test_field_srcs test/*.cpp)
    add_executable(test_field ${test_field_srcs})
    find_package(GTest CONFIG REQUIRED)
    target_link_libraries(test_field PRIVATE GTest::gtest GTest::gtest_main field)
    target_compile_definitions(
        test_field
        PRIVATE
        PEMU_MESH_TEST_DATA_DIR=\"${PROJECT_SOURCE_DIR}/submod/mesh/test\"
    )
    gtest_discover_tests(test_field)
endif()
