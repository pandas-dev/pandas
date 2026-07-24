#define PY_SSIZE_T_CLEAN
#include <stdint.h>
#include "CPy.h"

#define CPyFunction_weakreflist(f) (((PyCFunctionObject *)f)->m_weakreflist)
#define CPyFunction_class(f) ((PyObject*) ((PyCMethodObject *) (f))->mm_class)
#define CPyFunction_func_vectorcall(f) (((PyCFunctionObject *)f)->vectorcall)

static int
CPyFunction_clear(CPyFunction *m) {
    Py_CLEAR(((PyCFunctionObject*)m)->m_module);
    PyObject_ClearManagedDict((PyObject*)m);
    Py_CLEAR(m->func_name);
    Py_CLEAR(m->func_code);
    PyObject *cls = CPyFunction_class(m);
    ((PyCMethodObject *)m)->mm_class = NULL;
    Py_XDECREF(cls);

    return 0;
}

static void CPyFunction_dealloc(CPyFunction *m) {
    PyObject_GC_UnTrack(m);
    if (CPyFunction_weakreflist(m) != NULL)
        PyObject_ClearWeakRefs((PyObject *) m);
    CPyFunction_clear(m);
    PyMem_Free(m->func.func.m_ml);
    PyObject_GC_Del(m);
}

static PyObject* CPyFunction_repr(CPyFunction *op) {
    return PyUnicode_FromFormat("<function %U at %p>", op->func_name, (void *)op);
}

static PyObject* CPyFunction_call(PyObject *func, PyObject *args, PyObject *kw) {
    CPyFunction *f = (CPyFunction *)func;
    vectorcallfunc vc = CPyFunction_func_vectorcall(f);
    assert(vc);
    return PyVectorcall_Call(func, args, kw);
}

static int CPyFunction_traverse(CPyFunction *m, visitproc visit, void *arg) {
    Py_VISIT(((PyCFunctionObject *)m)->m_module);
    int e = PyObject_VisitManagedDict((PyObject*)m, visit, arg);
    if (e != 0) return e;
    Py_VISIT(m->func_name);
    Py_VISIT(m->func_code);
    Py_VISIT(CPyFunction_class(m));
    return 0;
}

static PyMemberDef CPyFunction_members[] = {
    {"__module__", T_OBJECT, offsetof(PyCFunctionObject, m_module), 0, 0},
    {"__vectorcalloffset__", T_PYSSIZET, offsetof(PyCFunctionObject, vectorcall), READONLY, 0},
    {"__weaklistoffset__", T_PYSSIZET, offsetof(PyCFunctionObject, m_weakreflist), READONLY, 0},
    {0, 0, 0,  0, 0}
};

PyObject* CPyFunction_get_name(PyObject *op, void *context) {
    (void)context;
    CPyFunction *func = (CPyFunction *)op;
    if (unlikely(func->func_name == NULL)) {
        func->func_name = PyUnicode_InternFromString(((PyCFunctionObject *)func)->m_ml->ml_name);
        if (unlikely(func->func_name == NULL))
            return NULL;
    }
    Py_INCREF(func->func_name);
    return func->func_name;
}

int CPyFunction_set_name(PyObject *op, PyObject *value, void *context) {
    (void)context;
    CPyFunction *func = (CPyFunction *)op;
    if (unlikely(!value || !PyUnicode_Check(value))) {
        PyErr_SetString(PyExc_TypeError, "__name__ must be set to a string object");
        return -1;
    }

    Py_INCREF(value);
    Py_XDECREF(func->func_name);
    func->func_name = value;
    return 0;
}

PyObject* CPyFunction_get_code(PyObject *op, void *context) {
    (void)context;
    CPyFunction *func = (CPyFunction *)op;
    PyObject* result = (func->func_code) ? func->func_code : Py_None;
    Py_INCREF(result);
    return result;
}

static PyObject* CPyFunction_get_none(PyObject *op, void *context) {
    (void)op;
    (void)context;
    PyObject* result = Py_None;
    Py_INCREF(result);
    return result;
}

int CPyFunction_set_none(PyObject *op, PyObject *value, void *context) {
    (void)op;
    (void)value;
    (void)context;
    return 0;
}

PyObject* CPyFunction_get_defaults(PyObject *op, void *context) {
    return CPyFunction_get_none(op, context);
}

PyObject* CPyFunction_get_kwdefaults(PyObject *op, void *context) {
    return CPyFunction_get_none(op, context);
}

PyObject* CPyFunction_get_annotations(PyObject *op, void *context) {
    return CPyFunction_get_none(op, context);
}

int CPyFunction_set_annotations(PyObject *op, PyObject *value, void *context) {
    return CPyFunction_set_none(op, value, context);
}

static PyGetSetDef CPyFunction_getsets[] = {
    {"__dict__", (getter)PyObject_GenericGetDict, (setter)PyObject_GenericSetDict, 0, 0},
    {"__name__", (getter)CPyFunction_get_name, (setter)CPyFunction_set_name, 0, 0},
    {"__code__", (getter)CPyFunction_get_code, 0, 0, 0},
    {"__defaults__", (getter)CPyFunction_get_defaults, 0, 0, 0},
    {"__kwdefaults__", (getter)CPyFunction_get_kwdefaults, 0, 0, 0},
    {"__annotations__", (getter)CPyFunction_get_annotations, CPyFunction_set_annotations, 0, 0},
    {0, 0, 0, 0, 0}
};

static PyObject* CPy_PyMethod_New(PyObject *func, PyObject *self, PyObject *typ) {
    (void)typ;
    if (!self) {
        Py_INCREF(func);
        return func;
    }
    return PyMethod_New(func, self);
}

static PyType_Slot CPyFunction_slots[] = {
    {Py_tp_dealloc, (void *)CPyFunction_dealloc},
    {Py_tp_repr, (void *)CPyFunction_repr},
    {Py_tp_call, (void *)CPyFunction_call},
    {Py_tp_traverse, (void *)CPyFunction_traverse},
    {Py_tp_clear, (void *)CPyFunction_clear},
    {Py_tp_members, (void *)CPyFunction_members},
    {Py_tp_getset, (void *)CPyFunction_getsets},
    {Py_tp_descr_get, (void *)CPy_PyMethod_New},
    {0, 0},
};

static PyType_Spec CPyFunction_spec = {
    .name = "Function compiled with mypyc",
    .basicsize = sizeof(CPyFunction),
    .itemsize = 0,
    .flags = Py_TPFLAGS_IMMUTABLETYPE |
#if PY_VERSION_HEX >= 0x030C0000
             Py_TPFLAGS_MANAGED_DICT |
#endif
             Py_TPFLAGS_HAVE_VECTORCALL | Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_BASETYPE,
    .slots = CPyFunction_slots,
};

static PyTypeObject *CPyFunctionType = NULL;

static PyObject* CPyFunction_Vectorcall(PyObject *func, PyObject *const *args, size_t nargsf, PyObject *kwnames) {
    CPyFunction *f = (CPyFunction *)func;
    Py_ssize_t nargs = PyVectorcall_NARGS(nargsf);
    PyObject *self;
    PyCFunction meth = ((PyCFunctionObject *)f)->m_ml->ml_meth;

    self = ((PyCFunctionObject *)f)->m_self;
    if (!self) {
        self = args[0];
        args += 1;
        nargs -= 1;
    }
    return ((_PyCFunctionFastWithKeywords)(void(*)(void))meth)(self, args, nargs, kwnames);
}


static CPyFunction* CPyFunction_Init(CPyFunction *op, PyMethodDef *ml, PyObject* name,
                                     PyObject *module, PyObject* code, bool set_self) {
    PyCFunctionObject *cf = (PyCFunctionObject *)op;
    CPyFunction_weakreflist(op) = NULL;
    cf->m_ml = ml;
    cf->m_self = set_self ? (PyObject *) op : NULL;

    Py_XINCREF(module);
    cf->m_module = module;

    Py_INCREF(name);
    op->func_name = name;

    ((PyCMethodObject *)op)->mm_class = NULL;

    Py_XINCREF(code);
    op->func_code = code;

    CPyFunction_func_vectorcall(op) = CPyFunction_Vectorcall;
    return op;
}

static PyObject* CPyCode_New(const char *filename, const char *funcname, int first_line, int flags) {
    PyCodeObject *code_obj = PyCode_NewEmpty(filename, funcname, first_line);
    if (unlikely(!code_obj)) {
        return NULL;
    }
    code_obj->co_flags = flags;
    return (PyObject *)code_obj;
}

static PyMethodDef* CPyMethodDef_New(const char *name, PyCFunction func, int flags, const char *doc) {
    PyMethodDef *method = (PyMethodDef *)PyMem_Malloc(sizeof(PyMethodDef));
    if (unlikely(!method)) {
        return NULL;
    }
    method->ml_name = name;
    method->ml_meth = func;
    method->ml_flags = flags;
    method->ml_doc = doc;
    return method;
}

PyObject* CPyFunction_New(PyObject *module, const char *filename, const char *funcname,
                          PyCFunction func, int func_flags, const char *func_doc,
                          int first_line, int code_flags, bool has_self_arg) {
    PyMethodDef *method = NULL;
    PyObject *code = NULL, *op = NULL;
    bool set_self = false;

    if (!CPyFunctionType) {
        CPyFunctionType = (PyTypeObject *)PyType_FromSpec(&CPyFunction_spec);
        if (unlikely(!CPyFunctionType)) {
            goto err;
        }
    }

    method = CPyMethodDef_New(funcname, func, func_flags, func_doc);
    if (unlikely(!method)) {
        goto err;
    }
    code = CPyCode_New(filename, funcname, first_line, code_flags);
    if (unlikely(!code)) {
        goto err;
    }

    // Set m_self inside the function wrapper only if the wrapped function has no self arg
    // to pass m_self as the self arg when the function is called.
    // When the function has a self arg, it will come in the args vector passed to the
    // vectorcall handler.
    set_self = !has_self_arg;
    op = (PyObject *)CPyFunction_Init(PyObject_GC_New(CPyFunction, CPyFunctionType),
                                      method, PyUnicode_FromString(funcname), module,
                                      code, set_self);
    if (unlikely(!op)) {
        goto err;
    }
    PyObject_GC_Track(op);
    return op;

err:
    CPyError_OutOfMemory();
    if (method) {
        PyMem_Free(method);
    }
    return NULL;
}
