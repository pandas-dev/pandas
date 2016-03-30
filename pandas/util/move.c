#include <Python.h>

#define COMPILING_IN_PY2 (PY_VERSION_HEX <= 0x03000000)

#if !COMPILING_IN_PY2
/* alias this because it is not aliased in Python 3 */
#define PyString_CheckExact PyBytes_CheckExact
#define PyString_AS_STRING PyBytes_AS_STRING
#define PyString_GET_SIZE PyBytes_GET_SIZE
#endif  /* !COMPILING_IN_PY2 */

#ifndef Py_TPFLAGS_HAVE_GETCHARBUFFER
#define Py_TPFLAGS_HAVE_GETCHARBUFFER 0
#endif

#ifndef Py_TPFLAGS_HAVE_NEWBUFFER
#define Py_TPFLAGS_HAVE_NEWBUFFER 0
#endif

PyObject *badmove;  /* bad move exception class */

typedef struct {
    PyObject_HEAD
    /* the bytes that own the buffer we are mutating */
    PyObject *invalid_bytes;
} stolenbufobject;

PyTypeObject stolenbuf_type;  /* forward declare type */

static void
stolenbuf_dealloc(stolenbufobject *self)
{
    Py_DECREF(self->invalid_bytes);
    PyObject_Del(self);
}

static int
stolenbuf_getbuffer(stolenbufobject *self, Py_buffer *view, int flags)
{
    return PyBuffer_FillInfo(view,
                             (PyObject*) self,
                             (void*) PyString_AS_STRING(self->invalid_bytes),
                             PyString_GET_SIZE(self->invalid_bytes),
                             0,  /* not readonly */
                             flags);
}

#if COMPILING_IN_PY2

static Py_ssize_t
stolenbuf_getreadwritebuf(stolenbufobject *self, Py_ssize_t segment, void **out)
{
    if (segment != 0) {
        PyErr_SetString(PyExc_SystemError,
                        "accessing non-existent string segment");
        return -1;
    }
    *out = PyString_AS_STRING(self->invalid_bytes);
    return PyString_GET_SIZE(self->invalid_bytes);
}

static Py_ssize_t
stolenbuf_getsegcount(stolenbufobject *self, Py_ssize_t *len)
{
    if (len) {
        *len = PyString_GET_SIZE(self->invalid_bytes);
    }
    return 1;
}

PyBufferProcs stolenbuf_as_buffer = {
    (readbufferproc) stolenbuf_getreadwritebuf,
    (writebufferproc) stolenbuf_getreadwritebuf,
    (segcountproc) stolenbuf_getsegcount,
    (charbufferproc) stolenbuf_getreadwritebuf,
    (getbufferproc) stolenbuf_getbuffer,
};

#else  /* Python 3 */

PyBufferProcs stolenbuf_as_buffer = {
    (getbufferproc) stolenbuf_getbuffer,
    NULL,
};

#endif  /* COMPILING_IN_PY2 */

static PyObject *
stolenbuf_new(PyObject *self, PyObject *args, PyObject *kwargs)
{
    stolenbufobject *ret;
    PyObject *bytes_rvalue;

    if (kwargs && PyDict_Size(kwargs)) {
        PyErr_SetString(PyExc_TypeError,
                        "stolenbuf does not accept keyword arguments");
        return NULL;
    }

    if (PyTuple_GET_SIZE(args) != 1) {
        PyErr_SetString(PyExc_TypeError,
                        "stolenbuf requires exactly 1 positional argument");
        return NULL;

    }

    /* pull out the single, positional argument */
    bytes_rvalue = PyTuple_GET_ITEM(args, 0);

    if (!PyString_CheckExact(bytes_rvalue)) {
        PyErr_SetString(PyExc_TypeError,
                        "stolenbuf can only steal from bytes objects");
        return NULL;
    }

    if (Py_REFCNT(bytes_rvalue) != 1) {
        /* there is a reference other than the caller's stack */
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

PyDoc_STRVAR(
    stolenbuf_doc,
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
    "If you want to use this function you are probably wrong.\n");

PyTypeObject stolenbuf_type = {
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
    0,                                          /* tp_traverse */
    0,                                          /* tp_clear */
    0,                                          /* tp_richcompare */
    0,                                          /* tp_weaklistoffset */
    0,                                          /* tp_iter */
    0,                                          /* tp_iternext */
    0,                                          /* tp_methods */
    0,                                          /* tp_members */
    0,                                          /* tp_getset */
    0,                                          /* tp_base */
    0,                                          /* tp_dict */
    0,                                          /* tp_descr_get */
    0,                                          /* tp_descr_set */
    0,                                          /* tp_dictoffset */
    0,                                          /* tp_init */
    0,                                          /* tp_alloc */
    (newfunc) stolenbuf_new,                    /* tp_new */
};

#define MODULE_NAME "pandas.util._move"

#if !COMPILING_IN_PY2
PyModuleDef _move_module = {
    PyModuleDef_HEAD_INIT,
    MODULE_NAME,
    NULL,
    -1,
};
#endif  /* !COMPILING_IN_PY2 */

PyDoc_STRVAR(
    badmove_doc,
    "Exception used to indicate that a move was attempted on a value with\n"
    "more than a single reference.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "data : any\n"
    "    The data which was passed to ``_move_into_mutable_buffer``.\n"
    "\n"
    "See Also\n"
    "--------\n"
    "pandas.util._move.move_into_mutable_buffer\n");

PyMODINIT_FUNC
#if !COMPILING_IN_PY2
#define ERROR_RETURN NULL
PyInit__move(void)
#else
#define ERROR_RETURN
init_move(void)
#endif  /* !COMPILING_IN_PY2 */
{
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

#if !COMPILING_IN_PY2
    if (!(m = PyModule_Create(&_move_module)))
#else
    if (!(m = Py_InitModule(MODULE_NAME, NULL)))
#endif  /* !COMPILING_IN_PY2 */
    {
        return ERROR_RETURN;
    }

    if (PyModule_AddObject(m,
                           "move_into_mutable_buffer",
                           (PyObject*) &stolenbuf_type)) {
        Py_DECREF(m);
        return ERROR_RETURN;
    }

    if (PyModule_AddObject(m, "BadMove", badmove)) {
        Py_DECREF(m);
        return ERROR_RETURN;
    }

#if !COMPILING_IN_PY2
    return m;
#endif  /* !COMPILING_IN_PY2 */
}
