#define PY_SSIZE_T_CLEAN
#include <Python.h>

static int series_getbuffer(PyObject *exporter, Py_buffer *view, int flags) {
    PyObject *arr = PyObject_CallMethod(exporter, "to_numpy", NULL);
    int result = PyObject_GetBuffer(arr, view, flags);
    Py_DECREF(arr);
    
    return result;
}

void series_releasebuffer(PyObject *exporter, Py_buffer *view) {
    return PyBuffer_Release(view);
}

static PyBufferProcs series_buffer_procs = {
    .bf_getbuffer = (getbufferproc)series_getbuffer,
    .bf_releasebuffer = (releasebufferproc)series_releasebuffer,
};

static PyTypeObject SeriesBufferMixin = {
    PyVarObject_HEAD_INIT(NULL, 0).tp_name = "SeriesBufferMixin",
    .tp_doc = "Mixin for buffering Series",
    .tp_basicsize = sizeof(PyObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = PyType_GenericNew,
    .tp_as_buffer = &series_buffer_procs,
};

static PyModuleDef seriesbuffermodule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "series_buffer",
    .m_doc = "Extension module for buffered Series.",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_series_buffer(void) {
    PyObject *m;
    if (PyType_Ready(&SeriesBufferMixin) < 0)
        return NULL;

    m = PyModule_Create(&seriesbuffermodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&SeriesBufferMixin);
    if (PyModule_AddObject(m, "SeriesBufferMixin",
                           (PyObject *)&SeriesBufferMixin) < 0) {
        Py_DECREF(&SeriesBufferMixin);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
