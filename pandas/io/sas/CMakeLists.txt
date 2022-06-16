add_library(_sas SHARED sas.c)
target_include_directories(_sas PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(_sas  ${Python_LIBRARIES})
