set(VTK_SHORT_VERSION 9.4)

vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO Kitware/VTK
    REF v9.4.2
    SHA512 28d293f326dac1c0476c70be0f5cafcfa4d5cbe015e7b9e57468a0632c3ad9c53820ce6d2beefe298bb7244cda3c25b43c6ef2b9ed2f1116bb462ee07539d867
    HEAD_REF master
)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DBUILD_TESTING=OFF
        -DVTK_BUILD_TESTING=OFF
        -DVTK_BUILD_EXAMPLES=OFF
        -DVTK_BUILD_ALL_MODULES=OFF
        -DVTK_ENABLE_REMOTE_MODULES=OFF
        -DVTK_FORBID_DOWNLOADS=ON
        -DVTK_GROUP_ENABLE_StandAlone=DONT_WANT
        -DVTK_GROUP_ENABLE_Rendering=NO
        -DVTK_GROUP_ENABLE_Views=NO
        -DVTK_MODULE_ENABLE_VTK_IOHDF=YES
        -DVTK_SMP_IMPLEMENTATION_TYPE=Sequential
        -DVTK_USE_EXTERNAL=OFF
        -DVTK_WRAP_JAVA=OFF
        -DVTK_WRAP_PYTHON=OFF
        -DVTK_WRAP_SERIALIZATION=OFF
)

vcpkg_cmake_install()
vcpkg_copy_pdbs()
vcpkg_cmake_config_fixup(CONFIG_PATH lib/cmake/vtk-${VTK_SHORT_VERSION})

set(VTK_TOOLS
    tokenize
    vtkEncodeString-${VTK_SHORT_VERSION}
    vtkHashSource-${VTK_SHORT_VERSION}
    vtkParseJava-${VTK_SHORT_VERSION}
    vtkParseOGLExt-${VTK_SHORT_VERSION}
    vtkProbeOpenGLVersion-${VTK_SHORT_VERSION}
    vtkTestOpenGLVersion-${VTK_SHORT_VERSION}
    vtkWrapHierarchy-${VTK_SHORT_VERSION}
    vtkWrapJava-${VTK_SHORT_VERSION}
    vtkWrapPython-${VTK_SHORT_VERSION}
    vtkWrapPythonInit-${VTK_SHORT_VERSION}
    vtkWrapSerDes-${VTK_SHORT_VERSION}
    vtkWrapTcl-${VTK_SHORT_VERSION}
    vtkWrapTclInit-${VTK_SHORT_VERSION}
)

foreach(tool IN LISTS VTK_TOOLS)
    set(release_tool
        "${CURRENT_PACKAGES_DIR}/bin/${tool}${VCPKG_TARGET_EXECUTABLE_SUFFIX}")
    if(EXISTS "${release_tool}")
        file(INSTALL
            "${release_tool}"
            DESTINATION "${CURRENT_PACKAGES_DIR}/tools/vtk"
            USE_SOURCE_PERMISSIONS)
        file(REMOVE "${release_tool}")
    endif()

    set(debug_tool
        "${CURRENT_PACKAGES_DIR}/debug/bin/${tool}${VCPKG_TARGET_EXECUTABLE_SUFFIX}")
    if(EXISTS "${debug_tool}")
        file(REMOVE "${debug_tool}")
    endif()
endforeach()

if(VCPKG_LIBRARY_LINKAGE STREQUAL "static")
    file(REMOVE_RECURSE
        "${CURRENT_PACKAGES_DIR}/bin"
        "${CURRENT_PACKAGES_DIR}/debug/bin")
endif()

file(REMOVE_RECURSE
    "${CURRENT_PACKAGES_DIR}/debug/include"
    "${CURRENT_PACKAGES_DIR}/debug/share")

if(EXISTS "${CURRENT_PACKAGES_DIR}/share/licenses")
    file(RENAME
        "${CURRENT_PACKAGES_DIR}/share/licenses"
        "${CURRENT_PACKAGES_DIR}/share/${PORT}/licenses")
endif()

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/Copyright.txt")
