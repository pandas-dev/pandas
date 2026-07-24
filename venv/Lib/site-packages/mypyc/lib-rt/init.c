#include <Python.h>
#include "CPy.h"
#include "static_data.c"

struct ExcDummyStruct _CPy_ExcDummyStruct = { PyObject_HEAD_INIT(NULL) };
PyObject *_CPy_ExcDummy = (PyObject *)&_CPy_ExcDummyStruct;

// System-wide empty tuple constant
PyObject * __mypyc_empty_tuple__ = NULL;

// Because its dynamic linker is more restricted than linux/OS X,
// Windows doesn't allow initializing globals with values from
// other dynamic libraries. This means we need to initialize
// things at load time.
void CPy_Init(void) {
    _CPy_ExcDummyStruct.ob_base.ob_type = &PyBaseObject_Type;

    // Initialize system-wide empty tuple constant
    if (__mypyc_empty_tuple__ == NULL) {
        __mypyc_empty_tuple__ = PyTuple_New(0);
        if (!__mypyc_empty_tuple__) {
            CPyError_OutOfMemory();
        }
    }
}
