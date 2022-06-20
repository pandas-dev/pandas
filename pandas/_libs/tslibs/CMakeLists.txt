set(BASIC_LIBRARIES base ccalendar dtypes nattype strptime timezones)
foreach(LIB ${BASIC_LIBRARIES})
  add_custom_command(OUTPUT ${LIB}.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                             ${LIB}.pyx)
  Python_add_library(${LIB} MODULE WITH_SOABI ${LIB}.c)
  target_include_directories(${LIB} PUBLIC ${Python_INCLUDE_DIRS}
                                           ${Python_NumPy_INCLUDE_DIRS})
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
  Python_add_library(${LIB} MODULE WITH_SOABI ${LIB}.c src/datetime/np_datetime.c)
  target_include_directories(
    ${LIB} PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS}
                  "src/datetime")
endforeach()

add_custom_command(OUTPUT np_datetime.c COMMAND ${Python_EXECUTABLE} -m cython
                                                -3 np_datetime.pyx)
Python_add_library(np_datetime MODULE WITH_SOABI np_datetime.c src/datetime/np_datetime.c
                               src/datetime/np_datetime_strings.c)
target_include_directories(np_datetime PUBLIC ${Python_INCLUDE_DIRS}
                                              ${Python_NumPy_INCLUDE_DIRS})

add_custom_command(OUTPUT parsing.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                            parsing.pyx)
Python_add_library(parsing MODULE WITH_SOABI parsing.c "../src/parser/tokenizer.c")
target_include_directories(
  parsing PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS}
                 "../src/klib")
