add_library(indexers SHARED indexers.c)
target_include_directories(indexers PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(indexers ${Python_LIBRARIES})

add_library(aggregations SHARED aggregations.cpp)
target_include_directories(aggregations PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(aggregations ${Python_LIBRARIES})
