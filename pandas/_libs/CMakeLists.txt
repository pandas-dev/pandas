set(BASIC_LIBRARIES arrays groupby hashing indexing internals reduction ops properties reshape)
foreach(LIB ${BASIC_LIBRARIES})
  add_library(${LIB} SHARED ${LIB}.c)
  target_include_directories(${LIB} PRIVATE ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
  target_link_libraries(${LIB} ${Python_LIBRARIES})
endforeach()

set(KLIB_REQUIRING_LIBRARIES algos hashtable join)
foreach(LIB ${KLIB_REQUIRING_LIBRARIES})
  add_library(${LIB} SHARED ${LIB}.c)
  target_include_directories(${LIB} PRIVATE ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} "src/klib")
  target_link_libraries(${LIB} ${Python_LIBRARIES})
endforeach()

# dispatch comes from ops_dispatch so can't use loop above
add_library(dispatch SHARED ops_dispatch.c)
target_include_directories(dispatch PRIVATE ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(dispatch ${Python_LIBRARIES})
