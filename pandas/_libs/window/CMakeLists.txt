add_custom_command(OUTPUT indexers.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                             indexers.pyx)
Python_add_library(indexers MODULE WITH_SOABI indexers.c)
target_include_directories(indexers PUBLIC ${Python_INCLUDE_DIRS}
                                           ${Python_NumPy_INCLUDE_DIRS})

add_custom_command(
  OUTPUT aggregations.cpp COMMAND ${Python_EXECUTABLE} -m cython -3 --cplus
                                  aggregations.pyx)
Python_add_library(aggregations MODULE WITH_SOABI aggregations.cpp)
target_include_directories(aggregations PUBLIC ${Python_INCLUDE_DIRS}
                                               ${Python_NumPy_INCLUDE_DIRS})
