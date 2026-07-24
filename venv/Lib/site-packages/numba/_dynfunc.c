/*
 * Definition of Environment and Closure objects.
 * This module is included by _dynfuncmod.c and by pycc-compiled modules.
 */

#include "_pymodule.h"

#include <string.h>


// if python version is 3.13
#if ((PY_MAJOR_VERSION == 3) && ((PY_MINOR_VERSION == 13) || PY_MINOR_VERSION == 14))
    #include "pythoncapi_compat.h"
    #define _Py_IsFinalizing Py_IsFinalizing
#endif
/* NOTE: EnvironmentObject and ClosureObject must be kept in sync with
 * the definitions in numba/targets/base.py (EnvBody and ClosureBody).
 */

/*
 * EnvironmentObject hosts data needed for execution of compiled functions.
 */
typedef struct {
    PyObject_HEAD
    PyObject *globals;
    /* Assorted "constants" that are needed at runtime to execute
       the compiled function.  This can include frozen closure variables,
       lifted loops, etc. */
    PyObject *consts;
} EnvironmentObject;


static PyMemberDef env_members[] = {
    {"globals", T_OBJECT, offsetof(EnvironmentObject, globals), READONLY, NULL},
    {"consts", T_OBJECT, offsetof(EnvironmentObject, consts), READONLY, NULL},
    {NULL}  /* Sentinel */
};

static int
env_traverse(EnvironmentObject *env, visitproc visit, void *arg)
{
    Py_VISIT(env->globals);
    Py_VISIT(env->consts);
    return 0;
}

static int
env_clear(EnvironmentObject *env)
{
    Py_CLEAR(env->globals);
    Py_CLEAR(env->consts);
    return 0;
}

static void
env_dealloc(EnvironmentObject *env)
{
    PyObject_GC_UnTrack((PyObject *) env);
    env_clear(env);
    Py_TYPE(env)->tp_free((PyObject *) env);
}

static EnvironmentObject *
env_new_empty(PyTypeObject* type)
{
    return (EnvironmentObject *) PyType_GenericNew(type, NULL, NULL);
}

static PyObject *
env_new(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
    PyObject *globals;
    EnvironmentObject *env;
    static char *kwlist[] = {"globals", 0};

    if (!PyArg_ParseTupleAndKeywords(
            args, kwds, "O!:function", kwlist,
            &PyDict_Type, &globals))
        return NULL;

    env = env_new_empty(type);
    if (env == NULL)
        return NULL;
    Py_INCREF(globals);
    env->globals = globals;
    env->consts = PyList_New(0);
    if (!env->consts) {
        Py_DECREF(env);
        return NULL;
    }
    return (PyObject *) env;
}


static PyTypeObject EnvironmentType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_dynfunc.Environment",                  /* tp_name */
    sizeof(EnvironmentObject),               /* tp_basicsize */
    0,                                       /* tp_itemsize */
    (destructor) env_dealloc,                /* tp_dealloc */
    0,                                       /* tp_vectorcall_offset */
    0,                                       /* tp_getattr*/
    0,                                       /* tp_setattr */
    0,                                       /* tp_as_async */
    0,                                       /* tp_repr */
    0,                                       /* tp_as_number */
    0,                                       /* tp_as_sequence */
    0,                                       /* tp_as_mapping */
    0,                                       /* tp_hash */
    0,                                       /* tp_call */
    0,                                       /* tp_str */
    0,                                       /* tp_getattro */
    0,                                       /* tp_setattro */
    0,                                       /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, /* tp_flags */
    0,                                       /* tp_doc */
    (traverseproc) env_traverse,             /* tp_traverse */
    (inquiry) env_clear,                     /* tp_clear */
    0,                                       /* tp_richcompare */
    0,                                       /* tp_weaklistoffset */
    0,                                       /* tp_iter */
    0,                                       /* tp_iternext */
    0,                                       /* tp_methods */
    env_members,                             /* tp_members */
    0,                                       /* tp_getset */
    0,                                       /* tp_base */
    0,                                       /* tp_dict */
    0,                                       /* tp_descr_get */
    0,                                       /* tp_descr_set */
    0,                                       /* tp_dictoffset */
    0,                                       /* tp_init */
    0,                                       /* tp_alloc */
    env_new,                                 /* tp_new */
    0,                                       /* tp_free */
    0,                                       /* tp_is_gc */
    0,                                       /* tp_bases */
    0,                                       /* tp_mro */
    0,                                       /* tp_cache */
    0,                                       /* tp_subclasses */
    0,                                       /* tp_weaklist */
    0,                                       /* tp_del */
    0,                                       /* tp_version_tag */
    0,                                       /* tp_finalize */
    0,                                       /* tp_vectorcall */
#if (PY_MAJOR_VERSION == 3) && (PY_MINOR_VERSION >= 12)
/* This was introduced first in 3.12
 * https://github.com/python/cpython/issues/91051
 */
    0,                                           /* tp_watched */
#endif
#if (PY_MAJOR_VERSION == 3) && (PY_MINOR_VERSION >= 13)
/* This was introduced in 3.13
 * https://github.com/python/cpython/pull/114900
 */
    0,                                           /* tp_versions_used */
#endif

/* WARNING: Do not remove this, only modify it! It is a version guard to
 * act as a reminder to update this struct on Python version update! */
#if (PY_MAJOR_VERSION == 3)
#if ! (NB_SUPPORTED_PYTHON_MINOR)
#error "Python minor version is not supported."
#endif
#else
#error "Python major version is not supported."
#endif
/* END WARNING*/
};

/* A closure object is created for each call to make_function(), and stored
   as the resulting PyCFunction object's "self" pointer.  It points to an
   EnvironmentObject which is constructed during compilation.  This allows
   for two things:
       - lifetime management of dependent data (e.g. lifted loop dispatchers)
       - access to the execution environment by the compiled function
         (for example the globals module)
   */

/* Closure is a variable-sized object for binary compatibility with
   Generator (see below). */
#define CLOSURE_HEAD          \
    PyObject_VAR_HEAD         \
    EnvironmentObject *env;

typedef struct {
    CLOSURE_HEAD
    /* The dynamically-filled method definition for the PyCFunction object
       using this closure. */
    PyMethodDef def;
    /* Arbitrary object to keep alive during the closure's lifetime.
       (put a tuple to put several objects alive).
       In practice, this helps keep the LLVM module and its generated
       code alive. */
    PyObject *keepalive;
    PyObject *weakreflist;
} ClosureObject;


static int
closure_traverse(ClosureObject *clo, visitproc visit, void *arg)
{
    Py_VISIT(clo->env);
    Py_VISIT(clo->keepalive);
    return 0;
}

static void
closure_dealloc(ClosureObject *clo)
{
    PyObject_GC_UnTrack((PyObject *) clo);
    if (clo->weakreflist != NULL)
        PyObject_ClearWeakRefs((PyObject *) clo);
    PyObject_Free((void *) clo->def.ml_name);
    PyObject_Free((void *) clo->def.ml_doc);
    Py_XDECREF(clo->env);
    Py_XDECREF(clo->keepalive);
    Py_TYPE(clo)->tp_free((PyObject *) clo);
}

static PyTypeObject ClosureType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_dynfunc._Closure",                     /* tp_name */
    sizeof(ClosureObject),                   /* tp_basicsize */
    0,                                       /* tp_itemsize */
    (destructor) closure_dealloc,            /* tp_dealloc */
    0,                                       /* tp_vectorcall_offset */
    0,                                       /* tp_getattr */
    0,                                       /* tp_setattr */
    0,                                       /* tp_as_async */
    0,                                       /* tp_repr */
    0,                                       /* tp_as_number */
    0,                                       /* tp_as_sequence */
    0,                                       /* tp_as_mapping */
    0,                                       /* tp_hash */
    0,                                       /* tp_call */
    0,                                       /* tp_str */
    0,                                       /* tp_getattro */
    0,                                       /* tp_setattro */
    0,                                       /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC, /* tp_flags */
    0,                                       /* tp_doc */
    (traverseproc) closure_traverse,         /* tp_traverse */
    0,                                       /* tp_clear */
    0,                                       /* tp_richcompare */
    offsetof(ClosureObject, weakreflist),    /* tp_weaklistoffset */
    0,                                       /* tp_iter */
    0,                                       /* tp_iternext */
    0,                                       /* tp_methods */
    0,                                       /* tp_members */
    0,                                       /* tp_getset */
    0,                                       /* tp_base */
    0,                                       /* tp_dict */
    0,                                       /* tp_descr_get */
    0,                                       /* tp_descr_set */
    0,                                       /* tp_dictoffset */
    0,                                       /* tp_init */
    0,                                       /* tp_alloc */
    0,                                       /* tp_new */
    0,                                       /* tp_free */
    0,                                       /* tp_is_gc */
    0,                                       /* tp_bases */
    0,                                       /* tp_mro */
    0,                                       /* tp_cache */
    0,                                       /* tp_subclasses */
    0,                                       /* tp_weaklist */
    0,                                       /* tp_del */
    0,                                       /* tp_version_tag */
    0,                                       /* tp_finalize */
    0,                                       /* tp_vectorcall */
#if (PY_MAJOR_VERSION == 3) && (PY_MINOR_VERSION >= 12)
/* This was introduced first in 3.12
 * https://github.com/python/cpython/issues/91051
 */
    0,                                           /* tp_watched */
#endif
#if (PY_MAJOR_VERSION == 3) && (PY_MINOR_VERSION >= 13)
/* This was introduced in 3.13
 * https://github.com/python/cpython/pull/114900
 */
    0,                                           /* tp_versions_used */
#endif

/* WARNING: Do not remove this, only modify it! It is a version guard to
 * act as a reminder to update this struct on Python version update! */
#if (PY_MAJOR_VERSION == 3)
#if ! (NB_SUPPORTED_PYTHON_MINOR)
#error "Python minor version is not supported."
#endif
#else
#error "Python major version is not supported."
#endif
/* END WARNING*/
};


/* Return an owned piece of character data duplicating a Python string
   object's value. */
static char *
dup_string(PyObject *strobj)
{
    const char *tmp = NULL;
    char *str;
    tmp = PyString_AsString(strobj);
    if (tmp == NULL)
        return NULL;
    /* Using PyObject_Malloc allows this memory to be tracked for
       leaks. */
    str = PyObject_Malloc(strlen(tmp) + 1);
    if (str == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    strcpy(str, tmp);
    return str;
}

/* Create and initialize a new Closure object */
static ClosureObject *
closure_new(PyObject *name, PyObject *doc, PyCFunction fnaddr,
            EnvironmentObject *env, PyObject *keepalive)
{
    ClosureObject *clo = (ClosureObject *) PyType_GenericAlloc(&ClosureType, 0);
    if (clo == NULL)
        return NULL;

    clo->def.ml_name = dup_string(name);
    if (!clo->def.ml_name) {
        Py_DECREF(clo);
        return NULL;
    }
    clo->def.ml_meth = fnaddr;
    clo->def.ml_flags = METH_VARARGS | METH_KEYWORDS;
    clo->def.ml_doc = dup_string(doc);
    if (!clo->def.ml_doc) {
        Py_DECREF(clo);
        return NULL;
    }
    Py_INCREF(env);
    clo->env = env;
    Py_XINCREF(keepalive);
    clo->keepalive = keepalive;
    return clo;
}

/* Create a new PyCFunction object wrapping a closure defined by
   the given arguments. */
static PyObject *
pycfunction_new(PyObject *module, PyObject *name, PyObject *doc,
                PyCFunction fnaddr, EnvironmentObject *env, PyObject *keepalive)
{
    PyObject *funcobj;
    PyObject *modname = NULL;
    ClosureObject *closure = NULL;

    closure = closure_new(name, doc, fnaddr, env, keepalive);
    if (closure == NULL) goto FAIL;

    modname = PyObject_GetAttrString(module, "__name__");
    if (modname == NULL) goto FAIL;

    funcobj = PyCFunction_NewEx(&closure->def, (PyObject *) closure, modname);
    Py_DECREF(closure);
    Py_DECREF(modname);

    return funcobj;

FAIL:
    Py_XDECREF(closure);
    Py_XDECREF(modname);
    return NULL;
}

/*
 * Python-facing wrapper for Numba-compiled generator.
 * Note the Environment's offset inside the struct is the same as in the
 * Closure object.  This is required to simplify generation of Python wrappers.
 */

typedef void (*gen_finalizer_t)(void *);

typedef struct {
    CLOSURE_HEAD
    PyCFunctionWithKeywords nextfunc;
    gen_finalizer_t finalizer;
    PyObject *weakreflist;
    union {
        double dummy;   /* Force alignment */
        char state[0];
    };
} GeneratorObject;

static int
generator_traverse(GeneratorObject *gen, visitproc visit, void *arg)
{
    /* XXX this doesn't traverse the state, which can own references to
       PyObjects */
    Py_VISIT(gen->env);
    return 0;
}

static int
generator_clear(GeneratorObject *gen)
{
    if (gen->finalizer != NULL) {
        gen->finalizer(gen->state);
        gen->finalizer = NULL;
    }
    Py_CLEAR(gen->env);
    gen->nextfunc = NULL;
    return 0;
}

static void
generator_dealloc(GeneratorObject *gen)
{
    PyObject_GC_UnTrack((PyObject *) gen);
    if (gen->weakreflist != NULL)
        PyObject_ClearWeakRefs((PyObject *) gen);
    /* XXX The finalizer may be called after the LLVM module has been
       destroyed (typically at interpreter shutdown) */
    if (!_Py_IsFinalizing())
        if (gen->finalizer != NULL)
            gen->finalizer(gen->state);
    Py_XDECREF(gen->env);
    Py_TYPE(gen)->tp_free((PyObject *) gen);
}

static PyObject *
generator_iternext(GeneratorObject *gen)
{
    PyObject *res, *args;
    if (gen->nextfunc == NULL) {
        PyErr_SetString(PyExc_RuntimeError,
                        "cannot call next() on finalized generator");
        return NULL;
    }
    args = PyTuple_Pack(1, (PyObject *) gen);
    if (args == NULL)
        return NULL;
    res = (*gen->nextfunc)((PyObject *) gen, args, NULL);
    Py_DECREF(args);
    return res;
}

static PyTypeObject GeneratorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_dynfunc._Generator",                    /* tp_name*/
    offsetof(GeneratorObject, state),         /* tp_basicsize*/
    1,                                        /* tp_itemsize*/
    (destructor) generator_dealloc,           /* tp_dealloc*/
    0,                                        /* tp_vectorcall_offset*/
    0,                                        /* tp_getattr*/
    0,                                        /* tp_setattr*/
    0,                                        /* tp_as_async*/
    0,                                        /* tp_repr*/
    0,                                        /* tp_as_number*/
    0,                                        /* tp_as_sequence*/
    0,                                        /* tp_as_mapping*/
    0,                                        /* tp_hash */
    0,                                        /* tp_call*/
    0,                                        /* tp_str*/
    0,                                        /* tp_getattro*/
    0,                                        /* tp_setattro*/
    0,                                        /* tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC
                       | Py_TPFLAGS_BASETYPE, /* tp_flags*/
    0,                                        /* tp_doc */
    (traverseproc) generator_traverse,        /* tp_traverse */
    (inquiry) generator_clear,                /* tp_clear */
    0,                                        /* tp_richcompare */
    offsetof(GeneratorObject, weakreflist),   /* tp_weaklistoffset */
    PyObject_SelfIter,                        /* tp_iter */
    (iternextfunc) generator_iternext,        /* tp_iternext */
    0,                                        /* tp_methods */
    0,                                        /* tp_members */
    0,                                        /* tp_getset */
    0,                                        /* tp_base */
    0,                                        /* tp_dict */
    0,                                        /* tp_descr_get */
    0,                                        /* tp_descr_set */
    0,                                        /* tp_dictoffset */
    0,                                        /* tp_init */
    0,                                        /* tp_alloc */
    0,                                        /* tp_new */
    0,                                        /* tp_free */
    0,                                        /* tp_is_gc */
    0,                                        /* tp_bases */
    0,                                        /* tp_mro */
    0,                                        /* tp_cache */
    0,                                        /* tp_subclasses */
    0,                                        /* tp_weaklist */
    0,                                        /* tp_del */
    0,                                        /* tp_version_tag */
    0,                                        /* tp_finalize */
    0,                                        /* tp_vectorcall */
#if (PY_MAJOR_VERSION == 3) && (PY_MINOR_VERSION >= 12)
/* This was introduced first in 3.12
 * https://github.com/python/cpython/issues/91051
 */
    0,                                           /* tp_watched */
#endif
#if (PY_MAJOR_VERSION == 3) && (PY_MINOR_VERSION >= 13)
/* This was introduced in 3.13
 * https://github.com/python/cpython/pull/114900
 */
    0,                                           /* tp_versions_used */
#endif

/* WARNING: Do not remove this, only modify it! It is a version guard to
 * act as a reminder to update this struct on Python version update! */
#if (PY_MAJOR_VERSION == 3)
#if ! (NB_SUPPORTED_PYTHON_MINOR)
#error "Python minor version is not supported."
#endif
#else
#error "Python major version is not supported."
#endif
/* END WARNING*/
};

/* Dynamically create a new generator object */
static PyObject *
Numba_make_generator(Py_ssize_t gen_state_size,
                     void *initial_state,
                     PyCFunctionWithKeywords nextfunc,
                     gen_finalizer_t finalizer,
                     EnvironmentObject *env)
{
    GeneratorObject *gen;
    gen = (GeneratorObject *) PyType_GenericAlloc(&GeneratorType, gen_state_size);
    if (gen == NULL)
        return NULL;
    memcpy(gen->state, initial_state, gen_state_size);
    gen->nextfunc = nextfunc;
    Py_XINCREF(env);
    gen->env = env;
    gen->finalizer = finalizer;
    return (PyObject *) gen;
}

/* Initialization subroutine for use by modules including this */
static int
init_dynfunc_module(PyObject *module)
{
    if (PyType_Ready(&ClosureType))
        return -1;
    if (PyType_Ready(&EnvironmentType))
        return -1;
    if (PyType_Ready(&GeneratorType))
        return -1;
    return 0;
}
