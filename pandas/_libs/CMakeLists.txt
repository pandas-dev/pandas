add_custom_command(
  OUTPUT algos_common_helper.pxi
         algos_take_helper.pxi
         hashtable_class_helper.pxi
         hashtable_func_helper.pxi
         index_class_helper.pxi
         intervaltree.pxi
         khash_for_primitive_helper.pxi
         sparse_op_helper.pxi
  COMMAND ${Python_EXECUTABLE} generate_templates.py)

add_custom_command(
  OUTPUT algos.c
  COMMAND ${Python_EXECUTABLE} -m cython -3 algos.pyx
  DEPENDS algos_common_helper.pxi algos_take_helper.pxi)
Python_add_library(algos MODULE WITH_SOABI algos.c)
target_include_directories(algos PUBLIC ${Python_INCLUDE_DIRS}
                                        ${Python_NumPy_INCLUDE_DIRS} "src/klib")

# There is a khash header file in src/klib and a cython generated
# one in _libs. Depending on the build timing one of the other could
# get picked up, though unclear why we need both? If we stick to the
# non-generated version we can remove any DEPENDS khash.h
add_custom_command(
  OUTPUT khash.h
  COMMAND ${Python_EXECUTABLE} -m cython -3 khash.pxd
  DEPENDS hkash_for_primitive_helper.pxi
)

add_custom_command(
  OUTPUT hashtable.c
  COMMAND ${Python_EXECUTABLE} -m cython -3 hashtable.pyx
  DEPENDS hashtable_class_helper.pxi hashtable_func_helper.pxi)

Python_add_library(hashtable MODULE WITH_SOABI hashtable.c)
target_include_directories(
  hashtable PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS}
                   "src/klib" DEPENDS khash.h)

add_custom_command(
  OUTPUT index.c
  COMMAND ${Python_EXECUTABLE} -m cython -3 index.pyx
  DEPENDS index_class_helper.pxi)
Python_add_library(index MODULE WITH_SOABI index.c)
target_include_directories(
  index PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib"
               "./tslibs")

add_custom_command(
  OUTPUT interval.c
  COMMAND ${Python_EXECUTABLE} -m cython -3 interval.pyx
  DEPENDS intervaltree.pxi)
Python_add_library(interval MODULE WITH_SOABI interval.c)
target_include_directories(
  interval PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS}
                  "src/klib" "./tslibs")

add_custom_command(
  OUTPUT sparse.c
  COMMAND ${Python_EXECUTABLE} -m cython -3 sparse.pyx
  DEPENDS sparse_op_helper.pxi)
Python_add_library(sparse MODULE WITH_SOABI sparse.c)
target_include_directories(
  sparse PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib")

# set(HEADERS arrays dtypes hashtable lib missing util) foreach(HEADER
# ${HEADERS}) add_custom_command(OUTPUT ${HEADER}.h COMMAND ${Python_EXECUTABLE}
# -m cython -3 ${HEADER}.pxd) endforeach()

set(BASIC_LIBRARIES
    arrays
    groupby
    hashing
    indexing
    internals
    reduction
    ops
    ops_dispatch
    properties
    reshape
    testing
    writers)
foreach(LIB ${BASIC_LIBRARIES})
  add_custom_command(OUTPUT ${LIB}.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                             ${LIB}.pyx)
  Python_add_library(${LIB} MODULE WITH_SOABI ${LIB}.c)
  target_include_directories(${LIB} PUBLIC ${Python_INCLUDE_DIRS}
                                           ${Python_NumPy_INCLUDE_DIRS})
endforeach()

add_subdirectory("tslibs")

add_custom_command(OUTPUT tslib.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                          tslib.pyx)
Python_add_library(tslib ${LIB} MODULE WITH_SOABI tslib.c tslibs/src/datetime/np_datetime.c)
target_include_directories(tslib PUBLIC ${Python_INCLUDE_DIRS}
                                        ${Python_NumPy_INCLUDE_DIRS} "./tslibs")

add_custom_command(OUTPUT missing.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                            missing.pyx)
Python_add_library(missing MODULE WITH_SOABI missing.c)
target_include_directories(
  missing PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "./tslibs")

add_custom_command(OUTPUT lib.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                        lib.pyx)
Python_add_library(lib MODULE WITH_SOABI lib.c src/parser/tokenizer.c)
target_include_directories(
  lib PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib"
             "./tslibs")

add_custom_command(OUTPUT join.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                         join.pyx DEPENDS khash_for_primitive_helper.pxi)
Python_add_library(join MODULE WITH_SOABI join.c)
target_include_directories(join PUBLIC ${Python_INCLUDE_DIRS}
                                       ${Python_NumPy_INCLUDE_DIRS} "src/klib")

add_custom_command(OUTPUT parsers.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                            parsers.pyx)
Python_add_library(parsers MODULE WITH_SOABI parsers.c src/parser/tokenizer.c src/parser/io.c)
target_include_directories(
  parsers PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib"
                 "src")

Python_add_library(
  ujson MODULE WITH_SOABI
  src/ujson/python/ujson.c
  src/ujson/python/objToJSON.c
  src/ujson/python/date_conversions.c
  src/ujson/python/JSONtoObj.c
  src/ujson/lib/ultrajsonenc.c
  src/ujson/lib/ultrajsondec.c
  tslibs/src/datetime/np_datetime.c
  tslibs/src/datetime/np_datetime_strings.c)
target_include_directories(
  ujson PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS}
               src/ujson/python src/ujson/lib src/datetime)

add_subdirectory("window")
