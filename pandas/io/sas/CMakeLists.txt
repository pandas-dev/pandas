add_custom_command(OUTPUT _sas.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                         _sas.pyx)
Python_add_library(_sas MODULE WITH_SOABI _sas.c)
target_include_directories(_sas PUBLIC ${Python_INCLUDE_DIRS}
                                      ${Python_NumPy_INCLUDE_DIRS})
