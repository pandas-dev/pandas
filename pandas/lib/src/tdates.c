#include <Python.h>
#include <datetime.h>
#include <numpy/arrayobject.h>
#include <stdio.h>

static PyObject* isAllDates(PyObject *self, PyObject *args) {
    PyDateTime_IMPORT;

    PyObject *input;
    PyArrayIterObject *iter;
    PyObject *obj;


    if (PyArg_ParseTuple(args, "O", &input)) {
        if (!PyArray_Check(input)) {
            PyErr_SetString(PyExc_RuntimeError, "Input was not ndarray!");
            return NULL;
        }

        long size = PyArray_SIZE(input);
    
        if (size == 0) {
            Py_RETURN_FALSE;	
        }
    
        iter = (PyArrayIterObject *) PyArray_IterNew(input);
    
        while (iter->index < iter->size) {
            obj = PyArray_GETITEM(input, (void *) iter->dataptr);
            if (!PyDateTime_Check(obj)) {
                Py_DECREF(obj);
                Py_DECREF(iter);
    
                Py_RETURN_FALSE;
            }
    
            Py_DECREF(obj);
            PyArray_ITER_NEXT(iter);
        }
    
        Py_DECREF(iter);
        Py_RETURN_TRUE;
    }	
    return NULL;
}

static PyMethodDef tdatesMethods[] =
{
    { "isAllDates", isAllDates, METH_VARARGS, NULL },
    {  NULL, NULL, 0, NULL},
};

PyMODINIT_FUNC inittdates(void)
{
    (void) Py_InitModule("tdates", tdatesMethods);
    import_array();
}
