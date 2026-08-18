include_guard(GLOBAL)

function(pemu_register_public_headers TARGET MODULE_NAME INCLUDE_DIR)
    set(PEMU_INCLUDE_ROOT
        "${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}"
    )
    set(MODULE_LINK
        "${PEMU_INCLUDE_ROOT}/${MODULE_NAME}"
    )
    file(MAKE_DIRECTORY
        "${PEMU_INCLUDE_ROOT}"
    )
    if(EXISTS "${MODULE_LINK}" OR IS_SYMLINK "${MODULE_LINK}")
        file(REMOVE "${MODULE_LINK}")
    endif()
    file(CREATE_LINK
        "${INCLUDE_DIR}"
        "${MODULE_LINK}"
        SYMBOLIC
    )
    target_include_directories(${TARGET}
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )
endfunction()