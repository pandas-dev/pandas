add_custom_command(OUTPUT sas.c COMMAND ${Python_EXECUTABLE} -m cython -3
                                        sas.pyx)
add_library(_sas SHARED sas.c)
target_include_directories(_sas PUBLIC ${Python_INCLUDE_DIRS}
                                      ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(_sas LINK_INTERFACE_LIBRARIES ${Python_LIBRARIES})
