# cython: profile=False
# distutils: language = c++
# cython: embedsignature = True

cdef extern from 'pandas/do_import_numpy.h':
    pass

cdef extern from 'pandas/numpy_interop.h' namespace 'pandas':
    void import_numpy()

cdef extern from 'pandas/init.h' namespace 'pandas':
    void libpandas_init()

# Ensure libpandas can use NumPy C Array API
import_numpy()
libpandas_init()
