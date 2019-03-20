#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <string.h>

#include "../../../src/inline_helper.h"
#include "../../../src/parser/tokenizer.h"

#if PY_MAJOR_VERSION >= 3
    #define PY_STRING_CHECK(string) (PyUnicode_Check(string))
#else
    #define PY_STRING_CHECK(string) \
        (PyString_Check(string) || PyUnicode_Check(string))
#endif

int PANDAS_INLINE convert_and_set_item(PyObject *item, Py_ssize_t index,
                                       PyArrayObject *result,
                                       int keep_trivial_numbers) {
    int needs_decref = 0, do_convert = 1;
    if (item == NULL) {
        return 0;
    }
    if (keep_trivial_numbers) {
        // don't convert an integer if it's zero,
        // don't convert a float if it's zero or NaN
#if PY_MAJOR_VERSION >= 3
        if (PyLong_Check(item)) {
            PyLongObject* v = (PyLongObject*)item;
            switch (Py_SIZE(v)) {
            case 0:
                do_convert = 0;
                break;
            case 1:  // fallthrough
            case -1:
                if (v->ob_digit[0] == 0) {
                    do_convert = 0;
                }
            }
#else
        if (PyInt_CheckExact(item)) {
            if (((PyIntObject*)item)->ob_ival == 0) do_convert = 0;
#endif
        } else if (PyFloat_Check(item)) {
            double v = PyFloat_AS_DOUBLE(item);
            if (v == 0.0 || v != v) {
                do_convert = 0;
            }
        }
    }

    if (do_convert) {
        if (!PY_STRING_CHECK(item)) {
            PyObject *str_item = PyObject_Str(item);
            if (str_item == NULL) {
                return 0;
            }
            item = str_item;
            needs_decref = 1;
        }
    }
    if (PyArray_SETITEM(result, PyArray_GETPTR1(result, index), item) != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot set resulting item");
        if (needs_decref) Py_DECREF(item);
        return 0;
    }
    if (needs_decref) Py_DECREF(item);
    return 1;
}

static int put_object_as_string(PyObject* list, Py_ssize_t idx,
                                PyObject* item) {
    if (!PY_STRING_CHECK(item)) {
        PyObject* str_item = PyObject_Str(item);
        if (str_item == NULL) {
            return 0;
        }
        Py_DECREF(item);
        item = str_item;
    }
    return (PyList_SetItem(list, idx, item) == 0) ? 1 : 0;
}

static PyObject* free_arrays(PyObject** arrays, Py_ssize_t size) {
    PyObject** item = arrays;
    Py_ssize_t i;
    for (i = 0; i < size; ++i, ++item) Py_DECREF(*item);
    free(arrays);
    return NULL;
}

static PyObject* concat_date_cols(PyObject *self, PyObject *args,
                                  PyObject *kwds) {
    PyObject *sequence = NULL;
    PyObject *py_keep_trivial_numbers = NULL;
    PyArrayObject *result = NULL;
    Py_ssize_t sequence_size = 0;
    int keep_trivial_numbers;
    char* kwlist[] = {"", "keep_trivial_numbers", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", kwlist,
                                     &sequence, &py_keep_trivial_numbers)) {
        return NULL;
    }
    if (!PySequence_Check(sequence)) {
        PyErr_SetString(PyExc_TypeError, "argument must be sequence");
        return NULL;
    }
    keep_trivial_numbers = (py_keep_trivial_numbers != NULL) ? \
            PyObject_IsTrue(py_keep_trivial_numbers) : 0;

    sequence_size = PySequence_Size(sequence);
    if (sequence_size == -1) {
        return NULL;
    } else if (sequence_size == 0) {
        npy_intp dims[1];
        dims[0] = 0;
        result = (PyArrayObject*)PyArray_ZEROS(1, dims, NPY_OBJECT, 0);
        return (PyObject*)result;
    } else if (sequence_size == 1) {
        PyObject* array = PySequence_GetItem(sequence, 0);
        Py_ssize_t array_size;
        if (array == NULL) {
            return NULL;
        }

        array_size = PySequence_Size(array);
        if (array_size == -1) {
            Py_DECREF(array);
            return NULL;
        }

        {
            npy_intp dims[1];
            dims[0] = array_size;
            result = (PyArrayObject*)PyArray_ZEROS(1, dims, NPY_OBJECT, 0);
            if (result == NULL) {
                Py_DECREF(array);
                return NULL;
            }
        }

        if (PyArray_CheckExact(array)) {
            PyArrayObject *ndarray = (PyArrayObject*)array;
            Py_ssize_t i;
            for (i = 0; i < array_size; ++i) {
                PyObject *item = PyArray_GETITEM(ndarray,
                                                 PyArray_GETPTR1(ndarray, i));
                if (!convert_and_set_item(item, i, result,
                                          keep_trivial_numbers)) {
                    Py_DECREF(result);
                    Py_DECREF(array);
                    Py_DECREF(item);
                    return NULL;
                }
                Py_DECREF(item);
            }
        } else {
            PyObject* fast_array = PySequence_Fast(array,
                    "elements of input sequence must be sequence");
            Py_ssize_t i;
            if (fast_array == NULL) {
                Py_DECREF(result);
                Py_DECREF(array);
                // PySequence_Fast set message, which in second argument
                return NULL;
            }

            for (i = 0; i < array_size; ++i) {
                PyObject* item = PySequence_Fast_GET_ITEM(fast_array, i);
                if (!convert_and_set_item(item, i, result,
                                          keep_trivial_numbers)) {
                    Py_DECREF(result);
                    Py_DECREF(array);
                    Py_DECREF(fast_array);
                    return NULL;
                }
            }
            Py_DECREF(fast_array);
        }
        Py_DECREF(array);
        return (PyObject*)result;
    } else {
        size_t mem_size = sizeof(PyObject*) * sequence_size;
        PyObject **arrays = (PyObject**) malloc(mem_size);
        PyObject *array = NULL;
        PyObject **parray = NULL;
        PyObject *fast_array = NULL;
        PyObject *separator = NULL;
        PyObject *item = NULL;
        PyObject *list_to_join = NULL;
        Py_ssize_t min_array_size = 0;
        int all_numpy = 1;
        Py_ssize_t i;
        for (i = 0; i < sequence_size; ++i) {
            array = PySequence_GetItem(sequence, i);
            if (array == NULL) {
                return free_arrays(arrays, i);
            }
            if (PyArray_CheckExact(array)) {
                if (PyArray_NDIM((PyArrayObject*)array) != 1) {
                    PyErr_SetString(PyExc_ValueError,
                                    "ndarrays must be 1-dimentional");
                    return free_arrays(arrays, i);
                }
            } else {
                all_numpy = 0;
            }
            arrays[i] = array;
        }

        parray = arrays;
        if (all_numpy) {
            Py_ssize_t i;
            for (i = 0; i < sequence_size; ++i, ++parray) {
                Py_ssize_t array_size = PyArray_SIZE((PyArrayObject*)(*parray));

                if (array_size < 0) {
                    return free_arrays(arrays, sequence_size);
                }

                if (array_size < min_array_size || min_array_size == 0) {
                    min_array_size = array_size;
                }
            }
        } else {
            Py_ssize_t i;
            for (i = 0; i < sequence_size; ++i, ++parray) {
                Py_ssize_t array_size;
                fast_array = PySequence_Fast(*parray,
                        "elements of input sequence must be sequence");
                array_size = (fast_array != NULL) ? \
                        PySequence_Fast_GET_SIZE(fast_array) : -1;

                if (array_size < 0) {
                    Py_XDECREF(fast_array);
                    return free_arrays(arrays, sequence_size);
                }
                Py_DECREF(*parray);
                arrays[i] = fast_array;

                if (array_size < min_array_size || min_array_size == 0) {
                    min_array_size = array_size;
                }
            }
        }

        {
            npy_intp dims[1];
            dims[0] = min_array_size;
            result = (PyArrayObject*)PyArray_ZEROS(1, dims, NPY_OBJECT, 0);
            if (result == NULL) {
                return free_arrays(arrays, sequence_size);
            }
        }

        separator = PyUnicode_FromFormat(" ");
        if (separator == NULL) {
            Py_DECREF(result);
            return free_arrays(arrays, sequence_size);
        }
        list_to_join = PyList_New(sequence_size);
        for (i = 0; i < min_array_size; ++i) {
            PyObject *result_string = NULL;
            parray = arrays;
            if (all_numpy) {
                Py_ssize_t j;
                for (j = 0; j < sequence_size; ++j, ++parray) {
                    PyArrayObject* arr = (PyArrayObject*)(*parray);
                    item = PyArray_GETITEM(arr, PyArray_GETPTR1(arr, i));
                    if (item == NULL) {
                        Py_DECREF(list_to_join);
                        Py_DECREF(result);
                        return free_arrays(arrays, sequence_size);
                    }
                    if (!put_object_as_string(list_to_join, j, item)) {
                        Py_DECREF(item);
                        Py_DECREF(list_to_join);
                        Py_DECREF(result);
                        return free_arrays(arrays, sequence_size);
                    }
                }
            } else {
                Py_ssize_t j;
                for (j = 0; j < sequence_size; ++j, ++parray) {
                    item = PySequence_Fast_GET_ITEM(*parray, i);
                    if (item == NULL) {
                        Py_DECREF(list_to_join);
                        Py_DECREF(result);
                        return free_arrays(arrays, sequence_size);
                    }
                    Py_INCREF(item);
                    if (!put_object_as_string(list_to_join, j, item)) {
                        Py_DECREF(item);
                        Py_DECREF(list_to_join);
                        Py_DECREF(result);
                        return free_arrays(arrays, sequence_size);
                    }
                }
            }
            result_string = PyUnicode_Join(separator, list_to_join);
            if (result_string == NULL) {
                Py_DECREF(list_to_join);
                Py_DECREF(result);
                return free_arrays(arrays, sequence_size);
            }
            if (PyArray_SETITEM(result, PyArray_GETPTR1(result, i),
                                result_string) != 0) {
                PyErr_SetString(PyExc_RuntimeError,
                                "Cannot set resulting item");
                Py_DECREF(list_to_join);
                Py_DECREF(result);
                Py_DECREF(result_string);
                return free_arrays(arrays, sequence_size);
            }
            Py_DECREF(result_string);
        }
        Py_DECREF(list_to_join);
        (void)free_arrays(arrays, sequence_size);
        return (PyObject*)result;
    }
}

static PyMethodDef module_methods[] = {
    /* name from python, name in C-file, ..., __doc__ string of method */
    {
        "concat_date_cols", (PyCFunction)concat_date_cols,
        METH_VARARGS | METH_KEYWORDS,
        "concatenates date cols and returns numpy array"
    },
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "datehelpers",                                   // name of module
    "helpers for datetime structures manipulation",  // module documentation
    -1,             // size of per-interpreter state of the module,
                    // or -1 if the module keeps state in global variables.
    module_methods
};
#define PY_DATEHELPERS_MODULE_INIT PyMODINIT_FUNC PyInit_datehelpers(void)
#define PY_MODULE_CREATE PyModule_Create(&moduledef)
#define PY_RETURN_MODULE return module
#else
#define PY_DATEHELPERS_MODULE_INIT void initdatehelpers(void)
#define PY_MODULE_CREATE Py_InitModule("datehelpers", module_methods)
#define PY_RETURN_MODULE
#endif

PY_DATEHELPERS_MODULE_INIT {
    PyObject *module = NULL;
    import_array();

    module = PY_MODULE_CREATE;

    PY_RETURN_MODULE;
}
