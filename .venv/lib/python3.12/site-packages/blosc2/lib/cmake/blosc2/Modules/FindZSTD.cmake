find_path(ZSTD_INCLUDE_DIR zstd.h)

find_library(ZSTD_LIBRARY NAMES zstd)

if(ZSTD_INCLUDE_DIR AND ZSTD_LIBRARY)
    set(ZSTD_FOUND TRUE)
    message(STATUS "Found ZSTD library: ${ZSTD_LIBRARY}")
endif()
