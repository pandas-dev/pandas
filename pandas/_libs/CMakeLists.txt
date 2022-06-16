set(BASIC_LIBRARIES arrays groupby hashing indexing internals reduction ops ops_dispatch properties reshape sparse testing writers)
foreach(LIB ${BASIC_LIBRARIES})
  add_library(${LIB} SHARED ${LIB}.c)
  target_include_directories(${LIB} PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
  target_link_libraries(${LIB} ${Python_LIBRARIES})
endforeach()

set(KLIB_REQUIRING algos hashtable join)
foreach(LIB ${KLIB_REQUIRING})
  add_library(${LIB} SHARED ${LIB}.c)
  target_include_directories(${LIB} PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib")
  target_link_libraries(${LIB} ${Python_LIBRARIES})
endforeach()

set(KLIB_AND_TSLIBS_REQUIRING index interval)
foreach(LIB ${KLIB_AND_TSLIBS_REQUIRING})
  add_library(${LIB} SHARED ${LIB}.c)
  target_include_directories(${LIB} PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib" "./tslibs")
  target_link_libraries(${LIB} ${Python_LIBRARIES})
endforeach()

SET(TSLIBS_REQUIRING missing tslib)
foreach(LIB ${TSLIBS_REQUIRING})
  add_library(${LIB} SHARED ${LIB}.c)
  target_include_directories(${LIB} PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "./tslibs")
  target_link_libraries(${LIB} ${Python_LIBRARIES})
endforeach()

add_subdirectory("tslibs")
add_subdirectory("window")

add_library(lib SHARED lib.c src/parser/tokenizer.c)
target_include_directories(lib PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib" "./tslibs")
target_link_libraries(lib ${Python_LIBRARIES})

add_library(parsers SHARED parsers.c src/parser/tokenizer.c src/parser/io.c)
target_include_directories(parsers PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib" "src")
target_link_libraries(parsers ${Python_LIBRARIES})

add_library(json SHARED src/ujson/python/ujson.c src/ujson/python/objToJSON.c
		 src/ujson/python/date_conversions.c src/ujson/python/JSONtoObj.c
		 src/ujson/lib/ultrajsonenc.c src/ujson/lib/ultrajsondec.c
		 tslibs/src/datetime/np_datetime.c tslibs/src/datetime/np_datetime_strings.c)
target_include_directories(json PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} src/ujson/python src/ujson/lib src/datetime)
target_link_libraries(json ${Python_LIBRARIES})
