file(GLOB_RECURSE output_srcs CONFIGURE_DEPENDS src/*.cpp)
add_library(output ${output_srcs})
add_library(pemu::output ALIAS output)
pemu_register_public_headers(
    output
    output
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

find_package(
    VTK 9.4 CONFIG REQUIRED
    COMPONENTS
        CommonCore
        CommonDataModel
        CommonExecutionModel
        FiltersCore
        IOHDF
)

target_link_libraries(
    output
    PUBLIC
        pemu::field
        pemu::mesh
        pemu::trace
        VTK::CommonCore
        VTK::CommonDataModel
        VTK::CommonExecutionModel
        VTK::IOHDF
)

vtk_module_autoinit(
    TARGETS output
    MODULES VTK::CommonCore VTK::CommonDataModel VTK::CommonExecutionModel VTK::IOHDF
)

if(BUILD_TESTING)
    file(GLOB_RECURSE test_output_srcs CONFIGURE_DEPENDS test/*.cpp)
    add_executable(test_output ${test_output_srcs})
    find_package(GTest CONFIG REQUIRED)
    target_link_libraries(
        test_output
        PRIVATE
            GTest::gtest
            GTest::gtest_main
            output
            VTK::IOHDF
    )
    target_compile_definitions(
        test_output
        PRIVATE
            PEMU_MESH_TEST_DATA_DIR="${PEMU_MESH_TEST_DATA_DIR}"
    )
    vtk_module_autoinit(TARGETS test_output MODULES VTK::IOHDF)
    gtest_discover_tests(test_output)
endif()
