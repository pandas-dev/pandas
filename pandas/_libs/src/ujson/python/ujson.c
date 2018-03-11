/*
Copyright (c) 2011-2013, ESN Social Software AB and Jonas Tarnstrom
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
* Neither the name of the ESN Social Software AB nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ESN SOCIAL SOFTWARE AB OR JONAS TARNSTROM BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Portions of code from MODP_ASCII - Ascii transformations (upper/lower, etc)
https://github.com/client9/stringencoders
Copyright (c) 2007  Nick Galbreath -- nickg [at] modp [dot] com. All rights reserved.

Numeric decoder derived from from TCL library
http://www.opensource.apple.com/source/tcl/tcl-14/tcl/license.terms
* Copyright (c) 1988-1993 The Regents of the University of California.
* Copyright (c) 1994 Sun Microsystems, Inc.
*/

#include "py_defines.h"
#include "version.h"

/* objToJSON */
PyObject *objToJSON(PyObject *self, PyObject *args, PyObject *kwargs);
void initObjToJSON(void);

/* JSONToObj */
PyObject *JSONToObj(PyObject *self, PyObject *args, PyObject *kwargs);

/* objToJSONFile */
PyObject *objToJSONFile(PyObject *self, PyObject *args, PyObject *kwargs);

/* JSONFileToObj */
PyObject *JSONFileToObj(PyObject *self, PyObject *args, PyObject *kwargs);

#define ENCODER_HELP_TEXT                                                  \
    "Use ensure_ascii=false to output UTF-8. Pass in double_precision to " \
    "alter the maximum digit precision of doubles. Set "                   \
    "encode_html_chars=True to encode < > & as unicode escape sequences."

static PyMethodDef ujsonMethods[] = {
    {"encode", (PyCFunction)objToJSON, METH_VARARGS | METH_KEYWORDS,
     "Converts arbitrary object recursively into JSON. " ENCODER_HELP_TEXT},
    {"decode", (PyCFunction)JSONToObj, METH_VARARGS | METH_KEYWORDS,
     "Converts JSON as string to dict object structure. Use precise_float=True "
     "to use high precision float decoder."},
    {"dumps", (PyCFunction)objToJSON, METH_VARARGS | METH_KEYWORDS,
     "Converts arbitrary object recursively into JSON. " ENCODER_HELP_TEXT},
    {"loads", (PyCFunction)JSONToObj, METH_VARARGS | METH_KEYWORDS,
     "Converts JSON as string to dict object structure. Use precise_float=True "
     "to use high precision float decoder."},
    {"dump", (PyCFunction)objToJSONFile, METH_VARARGS | METH_KEYWORDS,
     "Converts arbitrary object recursively into JSON "
     "file. " ENCODER_HELP_TEXT},
    {"load", (PyCFunction)JSONFileToObj, METH_VARARGS | METH_KEYWORDS,
     "Converts JSON as file to dict object structure. Use precise_float=True "
     "to use high precision float decoder."},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_libjson",
    0,            /* m_doc */
    -1,           /* m_size */
    ujsonMethods, /* m_methods */
    NULL,         /* m_reload */
    NULL,         /* m_traverse */
    NULL,         /* m_clear */
    NULL          /* m_free */
};

#define PYMODINITFUNC PyMODINIT_FUNC PyInit_json(void)
#define PYMODULE_CREATE() PyModule_Create(&moduledef)
#define MODINITERROR return NULL

#else

#define PYMODINITFUNC PyMODINIT_FUNC initjson(void)
#define PYMODULE_CREATE() Py_InitModule("json", ujsonMethods)
#define MODINITERROR return

#endif

PYMODINITFUNC {
    PyObject *module;
    PyObject *version_string;

    initObjToJSON();
    module = PYMODULE_CREATE();

    if (module == NULL) {
        MODINITERROR;
    }

    version_string = PyString_FromString(UJSON_VERSION);
    PyModule_AddObject(module, "__version__", version_string);

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
