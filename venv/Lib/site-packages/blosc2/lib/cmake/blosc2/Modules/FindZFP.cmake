# Find the ZFP compression library.
#
# Defines:
#   ZFP_FOUND
#   ZFP_INCLUDE_DIRS
#   ZFP_LIBRARIES
#   zfp            imported target, when found
#
# This module is intentionally small; prefer the upstream zfp CMake package
# when available via find_package(zfp CONFIG).

find_path(ZFP_INCLUDE_DIR
    NAMES zfp.h)

find_library(ZFP_LIBRARY
    NAMES zfp)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ZFP
    REQUIRED_VARS ZFP_LIBRARY ZFP_INCLUDE_DIR)

if(ZFP_FOUND)
    set(ZFP_INCLUDE_DIRS "${ZFP_INCLUDE_DIR}")
    set(ZFP_LIBRARIES "${ZFP_LIBRARY}")
    if(NOT TARGET zfp)
        add_library(zfp UNKNOWN IMPORTED)
        set_target_properties(zfp PROPERTIES
            IMPORTED_LOCATION "${ZFP_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${ZFP_INCLUDE_DIR}")
    endif()
endif()

mark_as_advanced(ZFP_INCLUDE_DIR ZFP_LIBRARY)
