//////////////////// LoadDataclassesModule.proto /////////////////////

static PyObject *__Pyx_Load_dataclasses_Module(void);

//////////////////// LoadDataclassesModule /////////////////////

static PyObject *__Pyx_Load_dataclasses_Module(void) {
    return PyImport_Import(PYIDENT("dataclasses"));
}

//////////////////// DataclassesCallHelper.proto ////////////////////////

static PyObject* __Pyx_DataclassesCallHelper(PyObject *callable, PyObject *kwds); /* proto */

//////////////////// DataclassesCallHelper ////////////////////////

// The signature of a few of the dataclasses module functions has
// been expanded over the years. Cython always passes the full set
// of arguments from the most recent version we know of, so needs
// to remove any arguments that don't exist on earlier versions.

static int __Pyx_DataclassesCallHelper_FilterToDict(PyObject *callable, PyObject *kwds, PyObject *new_kwds, PyObject *args_list, int is_kwonly) {
    Py_ssize_t size, i;
    size = PySequence_Size(args_list);
    if (unlikely(size < 0)) return -1;

    for (i=0; i<size; ++i) {
        PyObject *key, *value;
        int setitem_result;
        key = PySequence_GetItem(args_list, i);
        if (!key) return -1;

        if (PyUnicode_Check(key) && (
                PyUnicode_CompareWithASCIIString(key, "self") == 0 ||
                // namedtuple constructor in fallback code
                PyUnicode_CompareWithASCIIString(key, "_cls") == 0)) {
            Py_DECREF(key);
            continue;
        }

        value = PyDict_GetItem(kwds, key);
        if (!value) {
            if (is_kwonly) {
                Py_DECREF(key);
                continue;
            } else {
                // The most likely reason for this is that Cython
                // hasn't kept up to date with the Python dataclasses module.
                // To be nice to our users, try not to fail, but ask them
                // to report a bug so we can keep up to date.
                value = Py_None;
                if (PyErr_WarnFormat(
                        PyExc_RuntimeWarning, 1,
                        "Argument %S not passed to %R. This is likely a bug in Cython so please report it.",
                        key, callable) == -1) {
                    Py_DECREF(key);
                    return -1;
                }
            }
        }
        Py_INCREF(value);
        setitem_result = PyDict_SetItem(new_kwds, key, value);
        Py_DECREF(key);
        Py_DECREF(value);
        if (setitem_result == -1) return -1;
    }
    return 0;
}

static PyObject* __Pyx_DataclassesCallHelper(PyObject *callable, PyObject *kwds) {
    PyObject *new_kwds=NULL, *result=NULL;
    PyObject *inspect;
    PyObject *args_list=NULL, *kwonly_args_list=NULL, *getfullargspec_result=NULL;

    // Going via inspect to work out what arguments to pass is unlikely to be the
    // fastest thing ever. However, it is compatible, and only happens once
    // at module-import time.
    inspect = PyImport_ImportModule("inspect");
    if (!inspect) goto bad;
    getfullargspec_result = PyObject_CallMethodObjArgs(inspect, PYUNICODE("getfullargspec"), callable, NULL);
    Py_DECREF(inspect);
    if (!getfullargspec_result) goto bad;
    args_list = PyObject_GetAttrString(getfullargspec_result, "args");
    if (!args_list) goto bad;
    kwonly_args_list = PyObject_GetAttrString(getfullargspec_result, "kwonlyargs");
    if (!kwonly_args_list) goto bad;

    new_kwds = PyDict_New();
    if (!new_kwds) goto bad;

    // copy over only those arguments that are in the specification
    if (__Pyx_DataclassesCallHelper_FilterToDict(callable, kwds, new_kwds, args_list, 0) == -1) goto bad;
    if (__Pyx_DataclassesCallHelper_FilterToDict(callable, kwds, new_kwds, kwonly_args_list, 1) == -1) goto bad;
    result = PyObject_Call(callable, EMPTY(tuple), new_kwds);
bad:
    Py_XDECREF(getfullargspec_result);
    Py_XDECREF(args_list);
    Py_XDECREF(kwonly_args_list);
    Py_XDECREF(new_kwds);
    return result;
}
