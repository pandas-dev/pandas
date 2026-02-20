#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Blosc2::blosc2_shared" for configuration "Release"
set_property(TARGET Blosc2::blosc2_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Blosc2::blosc2_shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "/tmp/tmp2mxcauw4/wheel/platlib/blosc2/lib/libblosc2.so.2.23.1"
  IMPORTED_SONAME_RELEASE "libblosc2.so.7"
  )

list(APPEND _cmake_import_check_targets Blosc2::blosc2_shared )
list(APPEND _cmake_import_check_files_for_Blosc2::blosc2_shared "/tmp/tmp2mxcauw4/wheel/platlib/blosc2/lib/libblosc2.so.2.23.1" )

# Import target "Blosc2::blosc2_static" for configuration "Release"
set_property(TARGET Blosc2::blosc2_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Blosc2::blosc2_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "ASM;C"
  IMPORTED_LOCATION_RELEASE "/tmp/tmp2mxcauw4/wheel/platlib/blosc2/lib/libblosc2.a"
  )

list(APPEND _cmake_import_check_targets Blosc2::blosc2_static )
list(APPEND _cmake_import_check_files_for_Blosc2::blosc2_static "/tmp/tmp2mxcauw4/wheel/platlib/blosc2/lib/libblosc2.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
