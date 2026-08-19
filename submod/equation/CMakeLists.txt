file(GLOB_RECURSE equation_srcs src/*cpp)
add_library(equation)
target_sources(equation PRIVATE ${equation_srcs})
add_library(pemu::equation ALIAS equation)
pemu_register_public_headers(
    equation
    equation
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(
    equation
    PUBLIC
    pemu::discretization
)

if(BUILD_TESTING)
    file(GLOB_RECURSE test_equation_srcs test/*.cpp)
    add_executable(test_equation ${test_equation_srcs})
    find_package(GTest CONFIG REQUIRED)
    target_link_libraries(test_equation PRIVATE GTest::gtest GTest::gtest_main equation)
    target_compile_definitions(
        test_equation
        PRIVATE
        PEMU_MESH_TEST_DATA_DIR="${PEMU_MESH_TEST_DATA_DIR}"
    )
    gtest_discover_tests(test_equation)
endif()
