# file(GLOB_RECURSE physics_srcs src/*cpp)
add_library(physics INTERFACE)
# target_sources(physics PRIVATE ${physics_srcs})
add_library(pemu::physics ALIAS physics)
pemu_register_interface_headers(
    physics
    physics
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
