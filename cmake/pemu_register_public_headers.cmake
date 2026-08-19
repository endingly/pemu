include_guard(GLOBAL)

function(pemu_register_public_headers TARGET MODULE_NAME INCLUDE_DIR)
    set(MODULE_INCLUDE_ROOT
        "${CMAKE_BINARY_DIR}/include/${MODULE_NAME}")

    set(MODULE_LINK
        "${MODULE_INCLUDE_ROOT}/pemu/${MODULE_NAME}")

    file(MAKE_DIRECTORY
        "${MODULE_INCLUDE_ROOT}/pemu")

    file(CREATE_LINK
        "${INCLUDE_DIR}"
        "${MODULE_LINK}"
        SYMBOLIC)

    target_include_directories(${TARGET}
        PUBLIC
        $<BUILD_INTERFACE:${MODULE_INCLUDE_ROOT}>
        $<INSTALL_INTERFACE:include>)
endfunction()

function(pemu_register_interface_headers TARGET MODULE_NAME INCLUDE_DIR)
      set(MODULE_INCLUDE_ROOT
      "${CMAKE_BINARY_DIR}/include/${MODULE_NAME}")

  set(MODULE_LINK
      "${MODULE_INCLUDE_ROOT}/pemu/${MODULE_NAME}")

  file(MAKE_DIRECTORY
      "${MODULE_INCLUDE_ROOT}/pemu")

  file(CREATE_LINK
      "${INCLUDE_DIR}"
      "${MODULE_LINK}"
      SYMBOLIC)

  target_include_directories(${TARGET}
      INTERFACE
      $<BUILD_INTERFACE:${MODULE_INCLUDE_ROOT}>
      $<INSTALL_INTERFACE:include>)
endfunction()