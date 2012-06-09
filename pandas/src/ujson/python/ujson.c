#include "py_defines.h"
#include "version.h"

/* objToJSON */
PyObject* objToJSON(PyObject* self, PyObject *args, PyObject *kwargs);
void initObjToJSON(void);

/* JSONToObj */
PyObject* JSONToObj(PyObject* self, PyObject *args, PyObject *kwargs);

/* objToJSONFile */
PyObject* objToJSONFile(PyObject* self, PyObject *args, PyObject *kwargs);

/* JSONFileToObj */
PyObject* JSONFileToObj(PyObject* self, PyObject *args, PyObject *kwargs);


static PyMethodDef ujsonMethods[] = {
    {"encode", (PyCFunction) objToJSON, METH_VARARGS | METH_KEYWORDS, "Converts arbitrary object recursivly into JSON. Use ensure_ascii=false to output UTF-8. Pass in double_precision to alter the maximum digit precision with doubles"},
    {"decode", (PyCFunction) JSONToObj, METH_VARARGS | METH_KEYWORDS, "Converts JSON as string to dict object structure"},
    {"dumps", (PyCFunction) objToJSON, METH_VARARGS | METH_KEYWORDS,  "Converts arbitrary object recursivly into JSON. Use ensure_ascii=false to output UTF-8"},
    {"loads", (PyCFunction) JSONToObj, METH_VARARGS | METH_KEYWORDS,  "Converts JSON as string to dict object structure"},
    {"dump", (PyCFunction) objToJSONFile, METH_VARARGS | METH_KEYWORDS, "Converts arbitrary object recursively into JSON file. Use ensure_ascii=false to output UTF-8"},
    {"load", (PyCFunction) JSONFileToObj, METH_VARARGS | METH_KEYWORDS, "Converts JSON as file to dict object structure"},
    {NULL, NULL, 0, NULL}       /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_ujson",
    0,              /* m_doc */
    -1,             /* m_size */
    ujsonMethods,   /* m_methods */
    NULL,           /* m_reload */
    NULL,           /* m_traverse */
    NULL,           /* m_clear */
    NULL            /* m_free */
};

#define PYMODINITFUNC       PyObject *PyInit__ujson(void)
#define PYMODULE_CREATE()   PyModule_Create(&moduledef)
#define MODINITERROR        return NULL

#else

#define PYMODINITFUNC       PyMODINIT_FUNC init_ujson(void)
#define PYMODULE_CREATE()   Py_InitModule("_ujson", ujsonMethods)
#define MODINITERROR        return

#endif

PYMODINITFUNC
{
    PyObject *module;
    PyObject *version_string;

    initObjToJSON();   
    module = PYMODULE_CREATE();

    if (module == NULL)
    {
        MODINITERROR;
    }

    version_string = PyString_FromString (UJSON_VERSION);
    PyModule_AddObject (module, "__version__", version_string);

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
