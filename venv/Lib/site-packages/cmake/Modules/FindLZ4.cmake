find_path(LZ4_INCLUDE_DIR lz4.h)

find_library(LZ4_LIBRARY NAMES lz4 liblz4)

if(LZ4_INCLUDE_DIR AND LZ4_LIBRARY)
    set(LZ4_FOUND TRUE)
    message(STATUS "Found LZ4 library: ${LZ4_LIBRARY}")
else()
    message(STATUS "No LZ4 library found.  Using internal sources.")
endif()
