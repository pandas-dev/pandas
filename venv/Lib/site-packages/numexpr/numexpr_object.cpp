/*********************************************************************
  Numexpr - Fast numerical array expression evaluator for NumPy.

      License: MIT
      Author:  See AUTHORS.txt

  See LICENSE.txt for details about copyright and rights to use.
**********************************************************************/

#include "module.hpp"
#include <structmember.h>

#include "numexpr_config.hpp"
#include "interpreter.hpp"
#include "numexpr_object.hpp"

static int
size_from_char(char c)
{
    switch (c) {
        case 'b': return sizeof(char);
        case 'i': return sizeof(int);
        case 'l': return sizeof(long long);
        case 'f': return sizeof(float);
        case 'd': return sizeof(double);
        case 'c': return 2*sizeof(double);
        case 's': return 0;  /* strings are ok but size must be computed */
        default:
            PyErr_SetString(PyExc_TypeError, "signature value not in 'bilfdcs'");
            return -1;
    }
}

static void
NumExpr_dealloc(NumExprObject *self)
{
    Py_XDECREF(self->signature);
    Py_XDECREF(self->tempsig);
    Py_XDECREF(self->constsig);
    Py_XDECREF(self->fullsig);
    Py_XDECREF(self->program);
    Py_XDECREF(self->constants);
    Py_XDECREF(self->input_names);
    PyMem_Del(self->mem);
    PyMem_Del(self->rawmem);
    PyMem_Del(self->memsteps);
    PyMem_Del(self->memsizes);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
NumExpr_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    NumExprObject *self = (NumExprObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
#define INIT_WITH(name, object) \
        self->name = object; \
        if (!self->name) { \
            Py_DECREF(self); \
            return NULL; \
        }

        INIT_WITH(signature, PyBytes_FromString(""));
        INIT_WITH(tempsig, PyBytes_FromString(""));
        INIT_WITH(constsig, PyBytes_FromString(""));
        INIT_WITH(fullsig, PyBytes_FromString(""));
        INIT_WITH(program, PyBytes_FromString(""));
        INIT_WITH(constants, PyTuple_New(0));
        Py_INCREF(Py_None);
        self->input_names = Py_None;
        self->mem = NULL;
        self->rawmem = NULL;
        self->memsteps = NULL;
        self->memsizes = NULL;
        self->rawmemsize = 0;
        self->n_inputs = 0;
        self->n_constants = 0;
        self->n_temps = 0;
#undef INIT_WITH
    }
    return (PyObject *)self;
}

#define CHARP(s) ((char *)(s))

static int
NumExpr_init(NumExprObject *self, PyObject *args, PyObject *kwds)
{
    int i, j, mem_offset;
    int n_inputs, n_constants, n_temps;
    PyObject *signature = NULL, *tempsig = NULL, *constsig = NULL;
    PyObject *fullsig = NULL, *program = NULL, *constants = NULL;
    PyObject *input_names = NULL, *o_constants = NULL;
    int *itemsizes = NULL;
    char **mem = NULL, *rawmem = NULL;
    npy_intp *memsteps;
    npy_intp *memsizes;
    int rawmemsize;
    static char *kwlist[] = {CHARP("signature"), CHARP("tempsig"),
			     CHARP("program"),  CHARP("constants"),
			     CHARP("input_names"), NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "SSS|OO", kwlist,
                                     &signature,
                                     &tempsig,
                                     &program, &o_constants,
                                     &input_names)) {
        return -1;
    }

    n_inputs = (int)PyBytes_Size(signature);
    n_temps = (int)PyBytes_Size(tempsig);

    if (o_constants) {
        if (!PySequence_Check(o_constants) ) {
                PyErr_SetString(PyExc_TypeError, "constants must be a sequence");
                return -1;
        }
        n_constants = (int)PySequence_Length(o_constants);
        if (!(constants = PyTuple_New(n_constants)))
            return -1;
        if (!(constsig = PyBytes_FromStringAndSize(NULL, n_constants))) {
            Py_DECREF(constants);
            return -1;
        }
        if (!(itemsizes = PyMem_New(int, n_constants))) {
            Py_DECREF(constants);
            Py_DECREF(constsig);
            return -1;
        }
        for (i = 0; i < n_constants; i++) {
            PyObject *o;
            if (!(o = PySequence_GetItem(o_constants, i))) { /* new reference */
                Py_DECREF(constants);
                Py_DECREF(constsig);
                PyMem_Del(itemsizes);
                return -1;
            }
            PyTuple_SET_ITEM(constants, i, o); /* steals reference */
            if (PyBool_Check(o)) {
                PyBytes_AS_STRING(constsig)[i] = 'b';
                itemsizes[i] = size_from_char('b');
                continue;
            }

            if (PyArray_IsScalar(o, Int32)) {
                PyBytes_AS_STRING(constsig)[i] = 'i';
                itemsizes[i] = size_from_char('i');
                continue;
            }

            if (PyArray_IsScalar(o, Int64)) {
                PyBytes_AS_STRING(constsig)[i] = 'l';
                itemsizes[i] = size_from_char('l');
                continue;
            }
            /* The Float32 scalars are the only ones that should reach here */
            if (PyArray_IsScalar(o, Float32)) {
                PyBytes_AS_STRING(constsig)[i] = 'f';
                itemsizes[i] = size_from_char('f');
                continue;
            }
            if (PyFloat_Check(o)) {
                /* Python float constants are double precision by default */
                PyBytes_AS_STRING(constsig)[i] = 'd';
                itemsizes[i] = size_from_char('d');
                continue;
            }
            if (PyComplex_Check(o)) {
                PyBytes_AS_STRING(constsig)[i] = 'c';
                itemsizes[i] = size_from_char('c');
                continue;
            }
            if (PyBytes_Check(o)) {
                PyBytes_AS_STRING(constsig)[i] = 's';
                itemsizes[i] = (int)PyBytes_GET_SIZE(o);
                continue;
            }
            PyErr_SetString(PyExc_TypeError, "constants must be of type bool/int/long/float/double/complex/bytes");
            Py_DECREF(constsig);
            Py_DECREF(constants);
            PyMem_Del(itemsizes);
            return -1;
        }
    } else {
        n_constants = 0;
        if (!(constants = PyTuple_New(0)))
            return -1;
        if (!(constsig = PyBytes_FromString(""))) {
            Py_DECREF(constants);
            return -1;
        }
    }

    fullsig = PyBytes_FromFormat("%c%s%s%s", get_return_sig(program),
        PyBytes_AS_STRING(signature), PyBytes_AS_STRING(constsig),
        PyBytes_AS_STRING(tempsig));
    if (!fullsig) {
        Py_DECREF(constants);
        Py_DECREF(constsig);
        PyMem_Del(itemsizes);
        return -1;
    }

    if (!input_names) {
        input_names = Py_None;
    }

    /* Compute the size of registers. We leave temps out (will be
       malloc'ed later on). */
    rawmemsize = 0;
    for (i = 0; i < n_constants; i++)
        rawmemsize += itemsizes[i];
    rawmemsize *= BLOCK_SIZE1;

    mem = PyMem_New(char *, 1 + n_inputs + n_constants + n_temps);
    rawmem = PyMem_New(char, rawmemsize);
    memsteps = PyMem_New(npy_intp, 1 + n_inputs + n_constants + n_temps);
    memsizes = PyMem_New(npy_intp, 1 + n_inputs + n_constants + n_temps);
    if (!mem || !rawmem || !memsteps || !memsizes) {
        Py_DECREF(constants);
        Py_DECREF(constsig);
        Py_DECREF(fullsig);
        PyMem_Del(itemsizes);
        PyMem_Del(mem);
        PyMem_Del(rawmem);
        PyMem_Del(memsteps);
        PyMem_Del(memsizes);
        return -1;
    }
    /*
       0                                                  -> output
       [1, n_inputs+1)                                    -> inputs
       [n_inputs+1, n_inputs+n_consts+1)                  -> constants
       [n_inputs+n_consts+1, n_inputs+n_consts+n_temps+1) -> temps
    */
    /* Fill in 'mem' and 'rawmem' for constants */
    mem_offset = 0;
    for (i = 0; i < n_constants; i++) {
        char c = PyBytes_AS_STRING(constsig)[i];
        int size = itemsizes[i];
        mem[i+n_inputs+1] = rawmem + mem_offset;
        mem_offset += BLOCK_SIZE1 * size;
        memsteps[i+n_inputs+1] = memsizes[i+n_inputs+1] = size;
        /* fill in the constants */
        if (c == 'b') {
            char *bmem = (char*)mem[i+n_inputs+1];
            char value = (char)PyLong_AsLong(PyTuple_GET_ITEM(constants, i));
            for (j = 0; j < BLOCK_SIZE1; j++) {
                bmem[j] = value;
            }
        } else if (c == 'i') {
            int *imem = (int*)mem[i+n_inputs+1];
            int value = (int)PyLong_AsLong(PyTuple_GET_ITEM(constants, i));
            for (j = 0; j < BLOCK_SIZE1; j++) {
                imem[j] = value;
            }
        } else if (c == 'l') {
            long long *lmem = (long long*)mem[i+n_inputs+1];
            long long value = PyLong_AsLongLong(PyTuple_GET_ITEM(constants, i));
            for (j = 0; j < BLOCK_SIZE1; j++) {
                lmem[j] = value;
            }
        } else if (c == 'f') {
            /* In this particular case the constant is in a NumPy scalar
             and in a regular Python object */
            float *fmem = (float*)mem[i+n_inputs+1];
            float value = PyArrayScalar_VAL(PyTuple_GET_ITEM(constants, i),
                                            Float);
            for (j = 0; j < BLOCK_SIZE1; j++) {
                fmem[j] = value;
            }
        } else if (c == 'd') {
            double *dmem = (double*)mem[i+n_inputs+1];
            double value = PyFloat_AS_DOUBLE(PyTuple_GET_ITEM(constants, i));
            for (j = 0; j < BLOCK_SIZE1; j++) {
                dmem[j] = value;
            }
        } else if (c == 'c') {
            double *cmem = (double*)mem[i+n_inputs+1];
            Py_complex value = PyComplex_AsCComplex(PyTuple_GET_ITEM(constants, i));
            for (j = 0; j < 2*BLOCK_SIZE1; j+=2) {
                cmem[j] = value.real;
                cmem[j+1] = value.imag;
            }
        } else if (c == 's') {
            char *smem = (char*)mem[i+n_inputs+1];
            char *value = PyBytes_AS_STRING(PyTuple_GET_ITEM(constants, i));
            for (j = 0; j < size*BLOCK_SIZE1; j+=size) {
                memcpy(smem + j, value, size);
            }
        }
    }
    /* This is no longer needed since no unusual item sizes appear
       in temporaries (there are no string temporaries). */
    PyMem_Del(itemsizes);

    /* Fill in 'memsteps' and 'memsizes' for temps */
    for (i = 0; i < n_temps; i++) {
        char c = PyBytes_AS_STRING(tempsig)[i];
        int size = size_from_char(c);
        memsteps[i+n_inputs+n_constants+1] = size;
        memsizes[i+n_inputs+n_constants+1] = size;
    }
    /* See if any errors occured (e.g., in size_from_char) or if mem_offset is wrong */
    if (PyErr_Occurred() || mem_offset != rawmemsize) {
        if (mem_offset != rawmemsize) {
            PyErr_Format(PyExc_RuntimeError, "mem_offset does not match rawmemsize");
        }
        Py_DECREF(constants);
        Py_DECREF(constsig);
        Py_DECREF(fullsig);
        PyMem_Del(mem);
        PyMem_Del(rawmem);
        PyMem_Del(memsteps);
        PyMem_Del(memsizes);
        return -1;
    }


    #define REPLACE_OBJ(arg) \
    {PyObject *tmp = self->arg; \
     self->arg = arg; \
     Py_XDECREF(tmp);}
    #define INCREF_REPLACE_OBJ(arg) {Py_INCREF(arg); REPLACE_OBJ(arg);}
    #define REPLACE_MEM(arg) {PyMem_Del(self->arg); self->arg=arg;}

    INCREF_REPLACE_OBJ(signature);
    INCREF_REPLACE_OBJ(tempsig);
    REPLACE_OBJ(constsig);
    REPLACE_OBJ(fullsig);
    INCREF_REPLACE_OBJ(program);
    REPLACE_OBJ(constants);
    INCREF_REPLACE_OBJ(input_names);
    REPLACE_MEM(mem);
    REPLACE_MEM(rawmem);
    REPLACE_MEM(memsteps);
    REPLACE_MEM(memsizes);
    self->rawmemsize = rawmemsize;
    self->n_inputs = n_inputs;
    self->n_constants = n_constants;
    self->n_temps = n_temps;

    #undef REPLACE_OBJ
    #undef INCREF_REPLACE_OBJ
    #undef REPLACE_MEM

    return check_program(self);
}

static PyMethodDef NumExpr_methods[] = {
    {"run", (PyCFunction) NumExpr_run, METH_VARARGS|METH_KEYWORDS, NULL},
    {NULL, NULL}
};

static PyMemberDef NumExpr_members[] = {
    {CHARP("signature"), T_OBJECT_EX, offsetof(NumExprObject, signature), READONLY, NULL},
    {CHARP("constsig"), T_OBJECT_EX, offsetof(NumExprObject, constsig), READONLY, NULL},
    {CHARP("tempsig"), T_OBJECT_EX, offsetof(NumExprObject, tempsig), READONLY, NULL},
    {CHARP("fullsig"), T_OBJECT_EX, offsetof(NumExprObject, fullsig), READONLY, NULL},

    {CHARP("program"), T_OBJECT_EX, offsetof(NumExprObject, program), READONLY, NULL},
    {CHARP("constants"), T_OBJECT_EX, offsetof(NumExprObject, constants),
     READONLY, NULL},
    {CHARP("input_names"), T_OBJECT, offsetof(NumExprObject, input_names), 0, NULL},
    {NULL},
};

PyTypeObject NumExprType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "numexpr.NumExpr",         /*tp_name*/
    sizeof(NumExprObject),     /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)NumExpr_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    (ternaryfunc)NumExpr_run,  /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "NumExpr objects",         /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */
    NumExpr_methods,           /* tp_methods */
    NumExpr_members,           /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)NumExpr_init,    /* tp_init */
    0,                         /* tp_alloc */
    NumExpr_new,               /* tp_new */
};
