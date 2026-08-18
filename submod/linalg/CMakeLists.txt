file(GLOB_RECURSE linalg_srcs src/*cpp)
add_library(linalg ${linalg_srcs})
target_include_directories(linalg PUBLIC include)
find_package(Eigen3 CONFIG REQUIRED)
find_package(CHOLMOD CONFIG REQUIRED)
find_package(UMFPACK CONFIG REQUIRED)
target_link_libraries(
    linalg
    PUBLIC 
    Eigen3::Eigen
    SuiteSparse::CHOLMOD_static
    SuiteSparse::UMFPACK_static
)

if(BUILD_TESTING)
    file(GLOB_RECURSE test_linalg_srcs test/*.cpp)
    add_executable(test_linalg ${test_linalg_srcs})
    find_package(GTest CONFIG REQUIRED)
    target_link_libraries(test_linalg PRIVATE GTest::gtest GTest::gtest_main linalg)
    gtest_discover_tests(test_linalg)
endif()
