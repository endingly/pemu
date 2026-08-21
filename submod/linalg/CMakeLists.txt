file(GLOB_RECURSE linalg_srcs src/*cpp)
add_library(linalg)
target_sources(linalg PRIVATE ${linalg_srcs})
add_library(pemu::linalg ALIAS linalg)
pemu_register_public_headers(
    linalg
    linalg
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
find_package(Eigen3 CONFIG REQUIRED)
find_package(CHOLMOD CONFIG REQUIRED)
find_package(UMFPACK CONFIG REQUIRED)
target_link_libraries(
    linalg
    PUBLIC 
    Eigen3::Eigen
    SuiteSparse::CHOLMOD_static
    SuiteSparse::UMFPACK_static
    pemu::trace
)

if(BUILD_TESTING)
    file(GLOB_RECURSE test_linalg_srcs test/*.cpp)
    add_executable(test_linalg ${test_linalg_srcs})
    find_package(GTest CONFIG REQUIRED)
    target_link_libraries(test_linalg PRIVATE GTest::gtest GTest::gtest_main linalg)
    gtest_discover_tests(test_linalg)
endif()
