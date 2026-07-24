find_path(LZ4_INCLUDE_DIR lz4.h)

find_library(LZ4_LIBRARY NAMES lz4 liblz4)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LZ4 DEFAULT_MSG LZ4_LIBRARY LZ4_INCLUDE_DIR)

if(LZ4_FOUND)
    set(LZ4_INCLUDE_DIRS "${LZ4_INCLUDE_DIR}")
    set(LZ4_LIBRARIES "${LZ4_LIBRARY}")
    if(NOT TARGET LZ4::lz4)
        add_library(LZ4::lz4 UNKNOWN IMPORTED)
        set_target_properties(LZ4::lz4 PROPERTIES
            IMPORTED_LOCATION "${LZ4_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${LZ4_INCLUDE_DIR}")
    endif()
endif()

mark_as_advanced(LZ4_INCLUDE_DIR LZ4_LIBRARY)
