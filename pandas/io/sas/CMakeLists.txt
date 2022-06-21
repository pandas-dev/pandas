add_custom_command(OUTPUT _sas.c COMMAND ${Python3_EXECUTABLE} -m cython -3
                                         _sas.pyx)
python3_add_library(_sas MODULE WITH_SOABI _sas.c)
target_include_directories(_sas PUBLIC ${Python3_INCLUDE_DIRS}
                                       ${Python3_NumPy_INCLUDE_DIRS})
