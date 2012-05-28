#include <Python.h>

#if PY_MAJOR_VERSION >= 3

#define PyInt_Check             PyLong_Check
#define PyInt_AS_LONG           PyLong_AsLong
#define PyInt_FromLong          PyLong_FromLong

#define PyString_Check          PyBytes_Check
#define PyString_GET_SIZE       PyBytes_GET_SIZE
#define PyString_AS_STRING      PyBytes_AS_STRING

#define PyString_FromString     PyUnicode_FromString

#endif
