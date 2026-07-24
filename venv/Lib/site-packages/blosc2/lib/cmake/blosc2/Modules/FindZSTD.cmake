find_path(ZSTD_INCLUDE_DIR zstd.h)

find_library(ZSTD_LIBRARY NAMES zstd)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ZSTD DEFAULT_MSG ZSTD_LIBRARY ZSTD_INCLUDE_DIR)

if(ZSTD_FOUND)
    set(ZSTD_INCLUDE_DIRS "${ZSTD_INCLUDE_DIR}")
    set(ZSTD_LIBRARIES "${ZSTD_LIBRARY}")
    if(NOT TARGET zstd::libzstd)
        add_library(zstd::libzstd UNKNOWN IMPORTED)
        set_target_properties(zstd::libzstd PROPERTIES
            IMPORTED_LOCATION "${ZSTD_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${ZSTD_INCLUDE_DIR}")
    endif()
endif()

mark_as_advanced(ZSTD_INCLUDE_DIR ZSTD_LIBRARY)
