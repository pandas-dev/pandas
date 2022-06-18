add_custom_command(OUTPUT indexers.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                             indexers.pyx)
add_library(indexers SHARED indexers.c)
target_include_directories(indexers PUBLIC ${Python_INCLUDE_DIRS}
                                           ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(indexers PUBLIC ${Python_LIBRARIES})
target_link_options(indexers PUBLIC ${Python_LINK_OPTIONS})

add_custom_command(
  OUTPUT aggregations.cpp COMMAND ${Python_EXECUTABLE} -m cython -3 --cplus
                                  aggregations.pyx)
add_library(aggregations SHARED aggregations.cpp)
target_include_directories(aggregations PUBLIC ${Python_INCLUDE_DIRS}
                                               ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(aggregations PUBLIC ${Python_LIBRARIES})
target_link_options(aggregations PUBLIC ${Python_LINK_OPTIONS})
