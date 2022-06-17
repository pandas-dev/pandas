add_custom_command(OUTPUT khash.h khash_for_primitive_helper.pxi COMMAND python -m cythonize_templates khash header)

set(HEADERS arrays dtypes hashtable lib missing util)
foreach(HEADER ${HEADERS})
  add_custom_command(OUTPUT ${HEADER}.h COMMAND python -m cython -3 ${HEADER}.pxd)
endforeach()

set(BASIC_LIBRARIES arrays groupby hashing indexing internals reduction ops ops_dispatch properties reshape testing writers)
foreach(LIB ${BASIC_LIBRARIES})
  add_custom_command(OUTPUT ${LIB}.c COMMAND python -m cython -3 ${LIB}.pyx)
  add_library(${LIB} SHARED ${LIB}.c)
  target_include_directories(${LIB} PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
  target_link_libraries(${LIB} PUBLIC ${Python_LIBRARIES})
endforeach()

add_custom_command(OUTPUT sparse.c sparse_op_helper.pxi COMMAND python -m cythonize_templates sparse lib)
add_library(sparse SHARED sparse.c)
target_include_directories(sparse PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib")
target_link_libraries(sparse PUBLIC ${Python_LIBRARIES})

set(KLIB_REQUIRING join)
foreach(LIB ${KLIB_REQUIRING})
  add_custom_command(OUTPUT ${LIB}.c COMMAND python -m cython -3 ${LIB}.pyx DEPENDS khash.h)
  add_library(${LIB} SHARED ${LIB}.c)
  target_include_directories(${LIB} PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib")
  target_link_libraries(${LIB} PUBLIC ${Python_LIBRARIES})
endforeach()

add_custom_command(OUTPUT algos.c algos_common_helper.pxi algos_take_helper.pxi COMMAND python -m cythonize_templates algos lib DEPENDS khash.h)
add_library(algos SHARED algos.c)
target_include_directories(algos PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib")
target_link_libraries(algos PUBLIC ${Python_LIBRARIES})

add_custom_command(OUTPUT hashtable.c hashtable_class_helper.pxi hashtable_func_helper.pxi COMMAND python -m cythonize_templates hashtable lib DEPENDS khash.h)
add_library(hashtable SHARED hashtable.c)
target_include_directories(hashtable PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib")
target_link_libraries(hashtable PUBLIC ${Python_LIBRARIES})

add_custom_command(OUTPUT index.c index_class_helper.pxi COMMAND python -m cythonize_templates index lib DEPENDS khash.h)
add_library(index SHARED index.c)
target_include_directories(index PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib" "./tslibs")
target_link_libraries(index PUBLIC ${Python_LIBRARIES})

add_custom_command(OUTPUT interval.c intervaltree.pxi COMMAND python -m cythonize_templates interval lib DEPENDS khash.h)
add_library(interval SHARED interval.c)
target_include_directories(interval PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib" "./tslibs")
target_link_libraries(interval PUBLIC ${Python_LIBRARIES})

add_custom_command(OUTPUT tslib.c COMMAND python -m cython -3 tslib.pyx)
add_library(tslib ${LIB} SHARED tslib.c tslibs/src/datetime/np_datetime.c)
target_include_directories(tslib PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "./tslibs")
target_link_libraries(tslib PUBLIC ${Python_LIBRARIES})

add_custom_command(OUTPUT missing.c COMMAND python -m cython -3 missing.pyx)
add_library(missing SHARED missing.c)
target_include_directories(missing PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "./tslibs")
target_link_libraries(missing PUBLIC ${Python_LIBRARIES})

add_custom_command(OUTPUT lib.c COMMAND python -m cython -3 lib.pyx DEPENDS khash.h)
add_library(lib SHARED lib.c src/parser/tokenizer.c)
target_include_directories(lib PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib" "./tslibs")
target_link_libraries(lib PUBLIC ${Python_LIBRARIES})

add_custom_command(OUTPUT parsers.c COMMAND python -m cython -3 parsers.pyx DEPENDS khash.h)
add_library(parsers SHARED parsers.c src/parser/tokenizer.c src/parser/io.c)
target_include_directories(parsers PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib" "src")
target_link_libraries(parsers PUBLIC ${Python_LIBRARIES})

add_library(json SHARED src/ujson/python/ujson.c src/ujson/python/objToJSON.c
		 src/ujson/python/date_conversions.c src/ujson/python/JSONtoObj.c
		 src/ujson/lib/ultrajsonenc.c src/ujson/lib/ultrajsondec.c
		 tslibs/src/datetime/np_datetime.c tslibs/src/datetime/np_datetime_strings.c)
target_include_directories(json PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} src/ujson/python src/ujson/lib src/datetime)
target_link_libraries(json PUBLIC ${Python_LIBRARIES})

add_subdirectory("tslibs")
add_subdirectory("window")
