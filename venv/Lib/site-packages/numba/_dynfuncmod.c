#include "_dynfunc.c"

/* Python-facing function to dynamically create a new C function object */
static PyObject*
make_function(PyObject *self, PyObject *args)
{
    PyObject *module, *fname, *fdoc, *fnaddrobj;
    void *fnaddr;
    EnvironmentObject *env;
    PyObject *keepalive;

    if (!PyArg_ParseTuple(args, "OOOOO!|O",
            &module, &fname, &fdoc, &fnaddrobj, &EnvironmentType, &env,
            &keepalive)) {
        return NULL;
    }

    fnaddr = PyLong_AsVoidPtr(fnaddrobj);
    if (fnaddr == NULL && PyErr_Occurred())
        return NULL;

    return pycfunction_new(module, fname, fdoc, fnaddr, env, keepalive);
}

static PyMethodDef ext_methods[] = {
#define declmethod(func) { #func , ( PyCFunction )func , METH_VARARGS , NULL }
    declmethod(make_function),
    { NULL },
#undef declmethod
};


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

#define declmethod(func) _declpointer(#func, &Numba_##func)

#define declpointer(ptr) _declpointer(#ptr, &ptr)

    declmethod(make_generator);

#undef declmethod
    return dct;
error:
    Py_XDECREF(dct);
    return NULL;
}

MOD_INIT(_dynfunc) {
    PyObject *m, *impl_info;

    MOD_DEF(m, "_dynfunc", "No docs", ext_methods)
    if (m == NULL)
        return MOD_ERROR_VAL;

    if (init_dynfunc_module(m))
        return MOD_ERROR_VAL;

    impl_info = Py_BuildValue(
        "{snsnsn}",
        "offsetof_closure_body", offsetof(ClosureObject, env),
        "offsetof_env_body", offsetof(EnvironmentObject, globals),
        "offsetof_generator_state", offsetof(GeneratorObject, state)
        );
    if (impl_info == NULL)
        return MOD_ERROR_VAL;
    PyModule_AddObject(m, "_impl_info", impl_info);

    Py_INCREF(&ClosureType);
    PyModule_AddObject(m, "_Closure", (PyObject *) (&ClosureType));
    Py_INCREF(&EnvironmentType);
    PyModule_AddObject(m, "Environment", (PyObject *) (&EnvironmentType));
    Py_INCREF(&GeneratorType);
    PyModule_AddObject(m, "_Generator", (PyObject *) (&GeneratorType));

    PyModule_AddObject(m, "c_helpers", build_c_helpers_dict());

    return MOD_SUCCESS_VAL(m);
}
