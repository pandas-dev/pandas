#include <Python.h>
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
	{NULL, NULL, 0, NULL}		/* Sentinel */
};



PyMODINIT_FUNC
init_ujson(void)
{
	PyObject *module;
	PyObject *version_string;

	initObjToJSON();
	module = Py_InitModule("_ujson", ujsonMethods);

	version_string = PyString_FromString (UJSON_VERSION);
	PyModule_AddObject (module, "__version__", version_string);
}
