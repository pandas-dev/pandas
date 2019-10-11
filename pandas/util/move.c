/*
Copyright (c) 2019, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#ifndef Py_TPFLAGS_HAVE_GETCHARBUFFER
#define Py_TPFLAGS_HAVE_GETCHARBUFFER 0
#endif  // Py_TPFLAGS_HAVE_GETCHARBUFFER

#ifndef Py_TPFLAGS_HAVE_NEWBUFFER
#define Py_TPFLAGS_HAVE_NEWBUFFER 0
#endif  // Py_TPFLAGS_HAVE_NEWBUFFER

static PyObject *badmove;  /* bad move exception class */

typedef struct {
    PyObject_HEAD
    /* the bytes that own the buffer we are mutating */
    PyObject *invalid_bytes;
} stolenbufobject;

static PyTypeObject stolenbuf_type;  /* forward declare type */

static void
stolenbuf_dealloc(stolenbufobject *self) {
    Py_DECREF(self->invalid_bytes);
    PyObject_Del(self);
}

static int
stolenbuf_getbuffer(stolenbufobject *self, Py_buffer *view, int flags) {
    return PyBuffer_FillInfo(view,
                             (PyObject*) self,
                             (void*) PyBytes_AS_STRING(self->invalid_bytes),
                             PyBytes_GET_SIZE(self->invalid_bytes),
                             0,  /* not readonly */
                             flags);
}

static PyBufferProcs stolenbuf_as_buffer = {
    (getbufferproc) stolenbuf_getbuffer,
    NULL,
};

PyDoc_STRVAR(stolenbuf_doc,
             "A buffer that is wrapping a stolen bytes object's buffer.");

static PyTypeObject stolenbuf_type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "pandas.util._move.stolenbuf",              /* tp_name */
    sizeof(stolenbufobject),                    /* tp_basicsize */
    0,                                          /* tp_itemsize */
    (destructor) stolenbuf_dealloc,             /* tp_dealloc */
    0,                                          /* tp_print */
    0,                                          /* tp_getattr */
    0,                                          /* tp_setattr */
    0,                                          /* tp_reserved */
    0,                                          /* tp_repr */
    0,                                          /* tp_as_number */
    0,                                          /* tp_as_sequence */
    0,                                          /* tp_as_mapping */
    0,                                          /* tp_hash */
    0,                                          /* tp_call */
    0,                                          /* tp_str */
    0,                                          /* tp_getattro */
    0,                                          /* tp_setattro */
    &stolenbuf_as_buffer,                       /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
    Py_TPFLAGS_HAVE_NEWBUFFER |
    Py_TPFLAGS_HAVE_GETCHARBUFFER,              /* tp_flags */
    stolenbuf_doc,                              /* tp_doc */
};

PyDoc_STRVAR(
    move_into_mutable_buffer_doc,
    "Moves a bytes object that is about to be destroyed into a mutable buffer\n"
    "without copying the data.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "bytes_rvalue : bytes with 1 refcount.\n"
    "    The bytes object that you want to move into a mutable buffer. This\n"
    "    cannot be a named object. It must only have a single reference.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "buf : stolenbuf\n"
    "    An object that supports the buffer protocol which can give a mutable\n"
    "    view of the data that was previously owned by ``bytes_rvalue``.\n"
    "\n"
    "Raises\n"
    "------\n"
    "BadMove\n"
    "    Raised when a move is attempted on an object with more than one\n"
    "    reference.\n"
    "\n"
    "Notes\n"
    "-----\n"
    "If you want to use this function you are probably wrong.\n"
    "\n"
    "Warning: Do not call this function through *unpacking. This can\n"
    "potentially trick the reference checks which may allow you to get a\n"
    "mutable reference to a shared string!\n"
    "\n");

/* This is implemented as a standalone function instead of the ``tp_new`` of
   ``stolenbuf`` because we need to create a function using the METH_O flag
   to support Python 3.6. In python 3.6, PyCFunction calls from python code now
   count the reference owned by the argument tuple. This would cause the object
   to have 2 references if used with a direct call like: ``stolenbuf(a)``;
   however, if called through *unpacking like ``stolenbuf(*(a,))`` it would
   only have the one reference (the tuple). */
static PyObject*
move_into_mutable_buffer(PyObject *self, PyObject *bytes_rvalue) {
    stolenbufobject *ret;

    if (!PyBytes_CheckExact(bytes_rvalue)) {
        PyErr_SetString(PyExc_TypeError,
                        "stolenbuf can only steal from bytes objects");
        return NULL;
    }

    if (Py_REFCNT(bytes_rvalue) != 1) {
        // there is a reference other than the caller's stack
        PyErr_SetObject(badmove, bytes_rvalue);
        return NULL;
    }

    if (!(ret = PyObject_New(stolenbufobject, &stolenbuf_type))) {
        return NULL;
    }

    /* store the original bytes object in a field that is not
       exposed to python */
    Py_INCREF(bytes_rvalue);
    ret->invalid_bytes = bytes_rvalue;
    return (PyObject*) ret;
}

static PyMethodDef methods[] = {
    {"move_into_mutable_buffer",
     (PyCFunction) move_into_mutable_buffer,
     METH_O,
     move_into_mutable_buffer_doc},
    {NULL},
};

#define MODULE_NAME "pandas.util._move"

static PyModuleDef move_module = {
    PyModuleDef_HEAD_INIT,
    MODULE_NAME,
    NULL,
    -1,
    methods,
};

PyDoc_STRVAR(
    badmove_doc,
    "Exception used to indicate that a move was attempted on a value with\n"
    "more than a single reference.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "data : any\n"
    "    The data which was passed to ``move_into_mutable_buffer``.\n"
    "\n"
    "See Also\n"
    "--------\n"
    "pandas.util._move.move_into_mutable_buffer\n");

PyMODINIT_FUNC
#define ERROR_RETURN NULL
PyInit__move(void) {
    PyObject *m;

    if (!(badmove = PyErr_NewExceptionWithDoc("pandas.util._move.BadMove",
                                              badmove_doc,
                                              NULL,
                                              NULL))) {
        return ERROR_RETURN;
    }

    if (PyType_Ready(&stolenbuf_type)) {
        return ERROR_RETURN;
    }

    if (!(m = PyModule_Create(&move_module))) {
        return ERROR_RETURN;
    }

    if (PyModule_AddObject(m,
                           "stolenbuf",
                           (PyObject*) &stolenbuf_type)) {
        Py_DECREF(m);
        return ERROR_RETURN;
    }

    if (PyModule_AddObject(m, "BadMove", badmove)) {
        Py_DECREF(m);
        return ERROR_RETURN;
    }

    return m;
}
