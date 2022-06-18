set(BASIC_LIBRARIES base ccalendar dtypes nattype strptime timezones)
foreach(LIB ${BASIC_LIBRARIES})
  add_custom_command(OUTPUT ${LIB}.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                             ${LIB}.pyx)
  add_library(${LIB} SHARED ${LIB}.c)
  target_include_directories(${LIB} PUBLIC ${Python_INCLUDE_DIRS}
                                           ${Python_NumPy_INCLUDE_DIRS})
  target_link_libraries(${LIB} LINK_INTERFACE_LIBRARIES ${Python_LIBRARIES})
endforeach()

set(NP_DATETIME_REQUIRING
    conversion
    fields
    offsets
    period
    timedeltas
    timestamps
    tzconversion
    vectorized)
foreach(LIB ${NP_DATETIME_REQUIRING})
  add_custom_command(OUTPUT ${LIB}.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                             ${LIB}.pyx)
  add_library(${LIB} SHARED ${LIB}.c src/datetime/np_datetime.c)
  target_include_directories(
    ${LIB} PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS}
                  "src/datetime")
  target_link_libraries(${LIB} LINK_INTERFACE_LIBRARIES ${Python_LIBRARIES})
endforeach()

add_custom_command(OUTPUT np_datetime.c COMMAND ${Python_EXECUTABLE} -m cython
                                                -3 np_datetime.pyx)
add_library(np_datetime SHARED np_datetime.c src/datetime/np_datetime.c
                               src/datetime/np_datetime_strings.c)
target_include_directories(np_datetime PUBLIC ${Python_INCLUDE_DIRS}
                                              ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(np_datetime LINK_INTERFACE_LIBRARIES ${Python_LIBRARIES})

add_custom_command(OUTPUT parsing.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                            parsing.pyx)
add_library(parsing SHARED parsing.c "../src/parser/tokenizer.c")
target_include_directories(
  parsing PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS}
                 "../src/klib")
target_link_libraries(parsing LINK_INTERFACE_LIBRARIES ${Python_LIBRARIES})
