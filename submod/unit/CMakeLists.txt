# file(GLOB_RECURSE unit_srcs src/*.cpp)
add_library(unit INTERFACE)
add_library(pemu::unit ALIAS unit)

find_package(mp-units CONFIG REQUIRED)

pemu_register_interface_headers(
    unit
    unit
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(
    unit
    INTERFACE
    mp-units::mp-units
)
