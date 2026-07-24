/*
Expose all functions as pointers in a dedicated C extension.
*/
#include "cext/cext.h"
/* Import _pymodule.h first, for a recent _POSIX_C_SOURCE */
#include "_pymodule.h"

#include <math.h>
#ifdef _MSC_VER
    #define false 0
    #define true 1
    #define bool int
#else
    #include <stdbool.h>
#endif

/*
Include C-extension here
*/
#include "cext/cext.h"

/* Numba C helpers */
#include "_helperlib.c"

static PyObject *
build_c_helpers_dict(void)
{
    PyObject *dct = PyDict_New();
    if (dct == NULL)
        goto error;

#define _declpointer(name, value) do {                 \
    PyObject *o = PyLong_FromVoidPtr(value);           \
    if (o == NULL) goto error;                         \
    if (PyDict_SetItemString(dct, name, o)) {          \
        Py_DECREF(o);                                  \
        goto error;                                    \
    }                                                  \
    Py_DECREF(o);                                      \
} while (0)

#define declmethod(func) _declpointer(#func, &numba_##func)

#define declpointer(ptr) _declpointer(#ptr, &numba_##ptr)

    declmethod(fixed_fmod);
    declmethod(fixed_fmodf);
    declmethod(set_fnclex);

    declmethod(sdiv);
    declmethod(srem);
    declmethod(udiv);
    declmethod(urem);
    declmethod(frexp);
    declmethod(frexpf);
    declmethod(ldexp);
    declmethod(ldexpf);
    declmethod(cpow);
    declmethod(cpowf);
    declmethod(erf);
    declmethod(erff);
    declmethod(erfc);
    declmethod(erfcf);
    declmethod(gamma);
    declmethod(gammaf);
    declmethod(lgamma);
    declmethod(lgammaf);
    declmethod(nextafter);
    declmethod(nextafterf);
    declmethod(complex_adaptor);
    declmethod(adapt_ndarray);
    declmethod(ndarray_new);
    declmethod(extract_record_data);
    declmethod(get_buffer);
    declmethod(adapt_buffer);
    declmethod(release_buffer);
    declmethod(extract_np_datetime);
    declmethod(create_np_datetime);
    declmethod(extract_np_timedelta);
    declmethod(create_np_timedelta);
    declmethod(recreate_record);
    declmethod(fptoui);
    declmethod(fptouif);
    declmethod(gil_ensure);
    declmethod(gil_release);
    declmethod(fatal_error);
    declmethod(py_type);
    declmethod(unpack_slice);
    declmethod(do_raise);
    declmethod(unpickle);
    declmethod(runtime_build_excinfo_struct);
    declmethod(attempt_nocopy_reshape);
    declmethod(get_pyobject_private_data);
    declmethod(set_pyobject_private_data);
    declmethod(reset_pyobject_private_data);

    /* BLAS / LAPACK */
    declmethod(xxgemm);
    declmethod(xxgemv);
    declmethod(xxdot);
    declmethod(xxgetrf);
    declmethod(ez_xxgetri);
    declmethod(xxpotrf);
    declmethod(ez_rgeev);
    declmethod(ez_cgeev);
    declmethod(ez_xxxevd);
    declmethod(ez_gesdd);
    declmethod(ez_geqrf);
    declmethod(ez_xxgqr);
    declmethod(ez_gelsd);
    declmethod(xgesv);
    declmethod(xxnrm2);

    /* PRNG support */
    declmethod(get_py_random_state);
    declmethod(get_np_random_state);
    declmethod(get_internal_random_state);
    declmethod(rnd_shuffle);
    declmethod(rnd_init);
    declmethod(poisson_ptrs);

    /* Unicode string support */
    declmethod(extract_unicode);
    declmethod(gettyperecord);
    declmethod(get_PyUnicode_ExtendedCase);

    /* for gdb breakpoint */
    declmethod(gdb_breakpoint);

    /* for dictionary support */
    declmethod(test_dict);
    declmethod(dict_new_sized);
    declmethod(dict_set_method_table);
    declmethod(dict_free);
    declmethod(dict_length);
    declmethod(dict_lookup);
    declmethod(dict_insert);
    declmethod(dict_insert_ez);
    declmethod(dict_delitem);
    declmethod(dict_popitem);
    declmethod(dict_iter_sizeof);
    declmethod(dict_iter);
    declmethod(dict_iter_next);
    declmethod(dict_dump);

    /* for list support */
    declmethod(test_list);
    declmethod(list_new);
    declmethod(list_set_method_table);
    declmethod(list_free);
    declmethod(list_base_ptr);
    declmethod(list_size_address);
    declmethod(list_length);
    declmethod(list_allocated);
    declmethod(list_is_mutable);
    declmethod(list_set_is_mutable);
    declmethod(list_setitem);
    declmethod(list_getitem);
    declmethod(list_append);
    declmethod(list_delitem);
    declmethod(list_delete_slice);
    declmethod(list_iter_sizeof);
    declmethod(list_iter);
    declmethod(list_iter_next);

    /* for set support */
    declmethod(test_set);
    declmethod(set_new);
    declmethod(set_new_sized);
    declmethod(set_set_method_table);
    declmethod(set_free);
    declmethod(set_length);
    declmethod(set_add);
    declmethod(set_contains);
    declmethod(set_discard);
    declmethod(set_iter_sizeof);
    declmethod(set_iter);
    declmethod(set_iter_next);
    declmethod(set_dump);

#define MATH_UNARY(F, R, A) declmethod(F);
#define MATH_BINARY(F, R, A, B) declmethod(F);
    #include "mathnames.h"
#undef MATH_UNARY
#undef MATH_BINARY

#undef declmethod
    return dct;
error:
    Py_XDECREF(dct);
    return NULL;
}


/*
 * Helper to deal with flushing stdout
 */
PyAPI_FUNC(void) _numba_flush_stdout(void) ;

void
_numba_flush_stdout(void) {
  fflush(stdout);
}


static PyMethodDef ext_methods[] = {
    { "rnd_get_state", (PyCFunction) _numba_rnd_get_state, METH_O, NULL },
    { "rnd_get_py_state_ptr", (PyCFunction) _numba_rnd_get_py_state_ptr, METH_NOARGS, NULL },
    { "rnd_get_np_state_ptr", (PyCFunction) _numba_rnd_get_np_state_ptr, METH_NOARGS, NULL },
    { "rnd_seed", (PyCFunction) _numba_rnd_seed, METH_VARARGS, NULL },
    { "rnd_set_state", (PyCFunction) _numba_rnd_set_state, METH_VARARGS, NULL },
    { "rnd_shuffle", (PyCFunction) _numba_rnd_shuffle, METH_O, NULL },
    { "_import_cython_function", (PyCFunction) _numba_import_cython_function, METH_VARARGS, NULL },
    { NULL },
};

/*
 * These functions are exported by the module's DLL, to exercise ctypes / cffi
 * without relying on libc availability (see https://bugs.python.org/issue23606)
 */

PyAPI_FUNC(double) _numba_test_sin(double x);
PyAPI_FUNC(double) _numba_test_cos(double x);
PyAPI_FUNC(double) _numba_test_exp(double x);
PyAPI_FUNC(void) _numba_test_vsquare(int n, double *x, double *out);
PyAPI_FUNC(double) _numba_test_funcptr(double (*func)(double));
PyAPI_FUNC(bool) _numba_test_boolean(void);

double _numba_test_sin(double x)
{
    return sin(x);
}

double _numba_test_cos(double x)
{
    return cos(x);
}

double _numba_test_exp(double x)
{
    return exp(x);
}

void _numba_test_vsquare(int n, double *x, double *out)
{
    int i;
    for (i = 0; i < n; i++)
        out[i] = pow(x[i], 2.0);
}

void _numba_test_vcube(int n, double *x, double *out)
{
    int i;
    for (i = 0; i < n; i++)
        out[i] = pow(x[i], 3.0);
}

double _numba_test_funcptr(double (*func)(double))
{
    return func(1.5);
}

bool _numba_test_boolean()
{
    return true;
}

MOD_INIT(_helperlib) {
    PyObject *m;
    MOD_DEF(m, "_helperlib", "No docs", ext_methods)
    if (m == NULL)
        return MOD_ERROR_VAL;

    import_array();

    PyModule_AddObject(m, "c_helpers", build_c_helpers_dict());
    PyModule_AddIntConstant(m, "long_min", LONG_MIN);
    PyModule_AddIntConstant(m, "long_max", LONG_MAX);
    PyModule_AddIntConstant(m, "py_buffer_size", sizeof(Py_buffer));
    PyModule_AddIntConstant(m, "py_gil_state_size", sizeof(PyGILState_STATE));
    PyModule_AddIntConstant(m, "py_unicode_1byte_kind", PyUnicode_1BYTE_KIND);
    PyModule_AddIntConstant(m, "py_unicode_2byte_kind", PyUnicode_2BYTE_KIND);
    PyModule_AddIntConstant(m, "py_unicode_4byte_kind", PyUnicode_4BYTE_KIND);
#if (PY_MAJOR_VERSION == 3)
#if ((PY_MINOR_VERSION == 10) || (PY_MINOR_VERSION == 11))
    PyModule_AddIntConstant(m, "py_unicode_wchar_kind", PyUnicode_WCHAR_KIND);
#endif
#endif
    numba_rnd_ensure_global_init();

    return MOD_SUCCESS_VAL(m);
}
