vcpkg_from_bitbucket(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO fathomteam/moab
    REF ${VERSION}
    SHA512 1b3bed63fb1f1be2d7769ce70cd5b8381468b52bd4585002e1bb956f7830f7ae86227e80cf2ebb4f68a9df5a25056204e810c57c5ccb5fb6a65484bd6b325882
    HEAD_REF master
)

# vcpkg_extract_source_archive(
# SOURCE_PATH
# ARCHIVE "${ARCHIVE}"
# )
vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"

    OPTIONS

    # -----------------------------------------------------
    # pemu currently only needs MOAB core mesh functionality.
    # Keep the initial backend deliberately minimal.
    # -----------------------------------------------------
    -DENABLE_MPI=OFF
    -DENABLE_HDF5=OFF
    -DENABLE_NETCDF=OFF

    -DENABLE_METIS=OFF
    -DENABLE_PARMETIS=OFF
    -DENABLE_ZOLTAN=OFF

    -DENABLE_TEMPESTREMAP=OFF
    -DENABLE_PYMOAB=OFF

    -DENABLE_FORTRAN=OFF
    -DENABLE_BLASLAPACK=OFF

    -DENABLE_TESTING=OFF

    MAYBE_UNUSED_VARIABLES
    ENABLE_FORTRAN
    ENABLE_BLASLAPACK
    ENABLE_TESTING
)

vcpkg_cmake_install()

#
# MOAB 5.6 revamped its CMake packaging.
#
# If the installed config lands under lib/cmake/MOAB,
# normalize it into vcpkg's share/moab package location.
#
if(EXISTS "${CURRENT_PACKAGES_DIR}/lib/cmake/MOAB")
    vcpkg_cmake_config_fixup(
        PACKAGE_NAME moab
        CONFIG_PATH lib/cmake/MOAB
    )
elseif(EXISTS "${CURRENT_PACKAGES_DIR}/lib/cmake/moab")
    vcpkg_cmake_config_fixup(
        PACKAGE_NAME moab
        CONFIG_PATH lib/cmake/moab
    )
endif()

vcpkg_copy_pdbs()

#
# Headers are identical for Debug/Release.
#
file(
    REMOVE_RECURSE
    "${CURRENT_PACKAGES_DIR}/debug/include"
)

#
# Avoid duplicated package metadata in debug/.
#
file(
    REMOVE_RECURSE
    "${CURRENT_PACKAGES_DIR}/debug/share"
)

vcpkg_install_copyright(
    FILE_LIST
    "${SOURCE_PATH}/LICENSE"
)