# Bootstrap thirdparty dependencies
if ("$ENV{DISABLE_NATIVE_TOOLCHAIN}" STREQUAL "")
  # Enable toolchain variable if the environment is setup
  set(NATIVE_TOOLCHAIN ON)
  message(STATUS "Toolchain build.")

  # If toolchain is not set, pick a directory
  if ("$ENV{NATIVE_TOOLCHAIN}" STREQUAL "")
    set(ENV{NATIVE_TOOLCHAIN} "${CMAKE_CURRENT_SOURCE_DIR}/toolchain")
  endif()

  # Set the environment variables for dependent versions
  set(ENV{GCC_VERSION} "4.9.2")
  set(ENV{GFLAGS_VERSION} "2.0")
  set(ENV{GLOG_VERSION} "0.3.3-p1")
  set(ENV{GOOGLETEST_VERSION} "20151222")

  set(ENV{GFLAGS_PREFIX} "$ENV{NATIVE_TOOLCHAIN}/gflags-$ENV{GFLAGS_VERSION}")
  set(ENV{GLOG_PREFIX} "$ENV{NATIVE_TOOLCHAIN}/glog-$ENV{GLOG_VERSION}")
  set(ENV{GOOGLETEST_PREFIX}
    "$ENV{NATIVE_TOOLCHAIN}/googletest-$ENV{GOOGLETEST_VERSION}")

  # Setting SYSTEM_GCC will use the toolchain dependencies compiled with the original
  # host's compiler.
  if ("$ENV{SYSTEM_GCC}" STREQUAL "")
    set(GCC_ROOT $ENV{NATIVE_TOOLCHAIN}/gcc-$ENV{GCC_VERSION})
    set(CMAKE_C_COMPILER ${GCC_ROOT}/bin/gcc)
    set(CMAKE_CXX_COMPILER ${GCC_ROOT}/bin/g++)
  endif()

  # If the toolchain directory does not yet exists, we assume that the dependencies
  # should be downloaded. If the download script is not available fail the
  # configuration.
  if (NOT IS_DIRECTORY $ENV{NATIVE_TOOLCHAIN})
    set(BOOTSTRAP_CMD "${BUILD_SUPPORT_DIR}/bootstrap_toolchain.py")
    # Download and unpack the dependencies
    message(STATUS "Downloading and extracting dependencies.")
    execute_process(COMMAND ${BOOTSTRAP_CMD} RESULT_VARIABLE BOOTSTRAP_RESULT)
    if (${BOOTSTRAP_RESULT} EQUAL 0)
      message(STATUS "Toolchain bootstrap complete.")
    else()
      message(FATAL_ERROR "Toolchain bootstrap failed.")
    endif()
  else()
    message(STATUS "Native toolchain picked up at $ENV{NATIVE_TOOLCHAIN}")
  endif()
else()
  set(NATIVE_TOOLCHAIN OFF)
  message(STATUS "Native toolchain was explicitly disabled using DISABLE_NATIVE_TOOLCHAIN.")
  message(STATUS "Assuming system search path for dependencies.")
endif()
