add_custom_command(OUTPUT _sas.c COMMAND python -m cython -3 _sas.pyx)
add_library(_sas SHARED sas.c)
target_include_directories(_sas PUBLIC ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(_sas PUBLIC ${Python_LIBRARIES})
