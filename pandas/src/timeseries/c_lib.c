#include "c_lib.h"
#include "numpy/arrayobject.h"

// Numpy UFUNCS
static PyObject *NP_ADD, *NP_MULTIPLY, *NP_SUBTRACT, *NP_SQRT,
                *NP_GREATER, *NP_GREATER_EQUAL;

/*********************************************************
** Convenience wrappers for numpy UFUNCS                **
*********************************************************/
PyObject*
np_add(PyObject *left_val, PyObject *right_val) {

    PyObject *result;

    result = PyObject_CallFunction(
                         NP_ADD, "OO",
                         (PyArrayObject*)left_val,
                         right_val);
    return result;
}



PyObject*
np_subtract(PyObject *left_val, PyObject *right_val) {

    PyObject *result;

    result = PyObject_CallFunction(
                         NP_SUBTRACT, "OO",
                         (PyArrayObject*)left_val,
                         right_val);
    return result;
}



PyObject*
np_multiply(PyObject *left_val, PyObject *right_val) {

    PyObject *result;

    result = PyObject_CallFunction(
                         NP_MULTIPLY, "OO",
                         (PyArrayObject*)left_val,
                         right_val);
    return result;
}



PyObject*
np_sqrt(PyObject *val) {
    return PyObject_CallFunction(NP_SQRT, "(O)", val);
}



int np_greater(PyObject *left_val, PyObject *right_val) {

    PyObject *temp;
    int result;

    temp = PyObject_CallFunction(
                         NP_GREATER, "OO",
                         (PyArrayObject*)left_val,
                         right_val);

    result = (int)PyInt_AsLong(temp);
    Py_DECREF(temp);
    return result;
}



int np_greater_equal(PyObject *left_val, PyObject *right_val) {

    PyObject *temp;
    int result;

    temp = PyObject_CallFunction(
                         NP_GREATER_EQUAL, "OO",
                         (PyArrayObject*)left_val,
                         right_val);

    result = (int)PyInt_AsLong(temp);
    Py_DECREF(temp);
    return result;
}



char *str_uppercase(char *str) {
    if (str) {
        int i, len=strlen(str);
        char *result;
        if((result = PyArray_malloc((len + 1)*sizeof(char))) == NULL) {
            return (char *)PyErr_NoMemory();
        }
        strcpy(result, str);

        for (i=0;i<len;i++) {
            switch(result[i])
            {
                case 'a': { result[i]='A'; break; }
                case 'b': { result[i]='B'; break; }
                case 'c': { result[i]='C'; break; }
                case 'd': { result[i]='D'; break; }
                case 'e': { result[i]='E'; break; }
                case 'f': { result[i]='F'; break; }
                case 'g': { result[i]='G'; break; }
                case 'h': { result[i]='H'; break; }
                case 'i': { result[i]='I'; break; }
                case 'j': { result[i]='J'; break; }
                case 'k': { result[i]='K'; break; }
                case 'l': { result[i]='L'; break; }
                case 'm': { result[i]='M'; break; }
                case 'n': { result[i]='N'; break; }
                case 'o': { result[i]='O'; break; }
                case 'p': { result[i]='P'; break; }
                case 'q': { result[i]='Q'; break; }
                case 'r': { result[i]='R'; break; }
                case 's': { result[i]='S'; break; }
                case 't': { result[i]='T'; break; }
                case 'u': { result[i]='U'; break; }
                case 'v': { result[i]='V'; break; }
                case 'w': { result[i]='W'; break; }
                case 'x': { result[i]='X'; break; }
                case 'y': { result[i]='Y'; break; }
                case 'z': { result[i]='Z'; break; }
            }
        }

        return result;
    } else { return NULL; }
}



char *str_replace(const char *s, const char *old, const char *new) {
    char *ret;
    int i, count = 0;
    size_t newlen = strlen(new);
    size_t oldlen = strlen(old);

    for (i = 0; s[i] != '\0'; i++) {
        if (strstr(&s[i], old) == &s[i]) {
           count++;
           i += oldlen - 1;
        }
    }

    ret = PyArray_malloc(i + 1 + count * (newlen - oldlen));
    if (ret == NULL) {return (char *)PyErr_NoMemory();}

    i = 0;
    while (*s) {
        if (strstr(s, old) == s) {
            strcpy(&ret[i], new);
            i += newlen;
            s += oldlen;
        } else {
            ret[i++] = *s++;
        }
    }
    ret[i] = '\0';

    return ret;
}



PyObject *
set_callback(PyObject *args, PyObject **callback)
{
    PyObject *result = NULL;
    PyObject *temp;

    if (PyArg_ParseTuple(args, "O:set_callback", &temp)) {

        if (!PyCallable_Check(temp)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be callable");
            return NULL;
        }

        Py_XINCREF(temp);        // Add a reference to new callback
        Py_XDECREF(*callback);  // Dispose of previous callback
        *callback = temp;       // Remember new callback
        // Boilerplate to return "None"
        Py_INCREF(Py_None);
        result = Py_None;
    }
    return result;
}



void import_c_lib(PyObject *m) {
    PyObject *ops_dict;

    import_array();

    ops_dict = PyArray_GetNumericOps();

    NP_ADD = PyDict_GetItemString(ops_dict, "add");
    NP_MULTIPLY = PyDict_GetItemString(ops_dict, "multiply");
    NP_SUBTRACT = PyDict_GetItemString(ops_dict, "subtract");
    NP_SQRT = PyDict_GetItemString(ops_dict, "sqrt");
    NP_GREATER = PyDict_GetItemString(ops_dict, "greater");
    NP_GREATER_EQUAL = PyDict_GetItemString(ops_dict, "greater_equal");

    Py_INCREF(NP_ADD);
    Py_INCREF(NP_MULTIPLY);
    Py_INCREF(NP_SUBTRACT);
    Py_INCREF(NP_SQRT);
    Py_INCREF(NP_GREATER);
    Py_INCREF(NP_GREATER_EQUAL);

    Py_DECREF(ops_dict);
}
