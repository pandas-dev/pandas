# Merge several static libraries into a single output archive.
#
# Required variables:
#   OUTPUT      archive to create/replace
#   INPUTS      input archives separated with '|'
# Optional variables:
#   AR          archiver path
#   RANLIB      ranlib path
#   SYSTEM_NAME CMAKE_SYSTEM_NAME
#   MSVC        TRUE for MSVC-like toolchains

if(NOT DEFINED OUTPUT OR OUTPUT STREQUAL "")
    message(FATAL_ERROR "merge_static_libs.cmake requires OUTPUT")
endif()
if(NOT DEFINED INPUTS OR INPUTS STREQUAL "")
    message(FATAL_ERROR "merge_static_libs.cmake requires INPUTS")
endif()

string(REPLACE "|" ";" _inputs "${INPUTS}")
foreach(_input IN LISTS _inputs)
    if(NOT EXISTS "${_input}")
        message(FATAL_ERROR "Static library to merge does not exist: ${_input}")
    endif()
endforeach()

get_filename_component(_out_dir "${OUTPUT}" DIRECTORY)
get_filename_component(_out_name_we "${OUTPUT}" NAME_WE)
get_filename_component(_out_ext "${OUTPUT}" EXT)
# Some Windows archivers, including llvm-lib, reject output archive names that
# do not use the conventional .lib suffix.  Keep the original extension last.
set(_tmp "${_out_dir}/${_out_name_we}.merged${_out_ext}")
file(REMOVE "${_tmp}")

if(MSVC)
    if(NOT AR)
        message(FATAL_ERROR "MSVC static library merging requires AR/lib.exe")
    endif()
    execute_process(
        COMMAND "${AR}" /NOLOGO "/OUT:${_tmp}" ${_inputs}
        RESULT_VARIABLE _result)
elseif(SYSTEM_NAME STREQUAL "Darwin")
    find_program(_libtool NAMES libtool xcrun-libtool)
    if(NOT _libtool)
        execute_process(
            COMMAND xcrun -find libtool
            OUTPUT_VARIABLE _libtool
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
    endif()
    if(NOT _libtool)
        message(FATAL_ERROR "Could not find libtool for static library merging")
    endif()
    execute_process(
        COMMAND "${_libtool}" -static -o "${_tmp}" ${_inputs}
        RESULT_VARIABLE _result)
else()
    if(NOT AR)
        message(FATAL_ERROR "Static library merging requires AR")
    endif()
    set(_mri "${_out_dir}/${_out_name_we}.mri")
    file(WRITE "${_mri}" "CREATE ${_tmp}\n")
    foreach(_input IN LISTS _inputs)
        file(APPEND "${_mri}" "ADDLIB ${_input}\n")
    endforeach()
    file(APPEND "${_mri}" "SAVE\nEND\n")
    execute_process(
        COMMAND "${AR}" -M
        INPUT_FILE "${_mri}"
        RESULT_VARIABLE _result)
    file(REMOVE "${_mri}")
endif()

if(NOT _result EQUAL 0)
    file(REMOVE "${_tmp}")
    message(FATAL_ERROR "Failed to merge static libraries into ${OUTPUT}")
endif()

if(RANLIB AND NOT MSVC)
    execute_process(COMMAND "${RANLIB}" "${_tmp}" RESULT_VARIABLE _ranlib_result)
    if(NOT _ranlib_result EQUAL 0)
        file(REMOVE "${_tmp}")
        message(FATAL_ERROR "ranlib failed for merged static library ${_tmp}")
    endif()
endif()

file(REMOVE "${OUTPUT}")
file(RENAME "${_tmp}" "${OUTPUT}")
