//////////////////// CythonFunctionShared.module_state_decls ////////////////////

PyTypeObject *__pyx_CyFunctionType;

//////////////////// CythonFunctionShared.module_state_traverse ///////////////////

Py_VISIT(traverse_module_state->__pyx_CyFunctionType);

//////////////////// CythonFunctionShared.module_state_clear ///////////////////

Py_CLEAR(clear_module_state->__pyx_CyFunctionType);

//////////////////// CythonFunctionShared.init //////////////////
//@substitute: naming

if (likely(__pyx_CyFunction_init($module_cname) == 0)); else

//////////////////// CythonFunctionShared.proto ////////////////////

#define __Pyx_CyFunction_USED

#define __Pyx_CYFUNCTION_STATICMETHOD  0x01
#define __Pyx_CYFUNCTION_CLASSMETHOD   0x02
#define __Pyx_CYFUNCTION_CCLASS        0x04
#define __Pyx_CYFUNCTION_COROUTINE     0x08

#define __Pyx_CyFunction_GetClosure(f) \
    (((__pyx_CyFunctionObject *) (f))->func_closure)

#if PY_VERSION_HEX < 0x030900B1 || CYTHON_COMPILING_IN_LIMITED_API
  #define __Pyx_CyFunction_GetClassObj(f) \
      (((__pyx_CyFunctionObject *) (f))->func_classobj)
#else
  #define __Pyx_CyFunction_GetClassObj(f) \
      ((PyObject*) ((PyCMethodObject *) (f))->mm_class)
#endif
#define __Pyx_CyFunction_SetClassObj(f, classobj)  \
    __Pyx__CyFunction_SetClassObj((__pyx_CyFunctionObject *) (f), (classobj))

#define __Pyx_CyFunction_Defaults(type, f) \
    ((type *)(((__pyx_CyFunctionObject *) (f))->defaults))
#define __Pyx_CyFunction_SetDefaultsGetter(f, g) \
    ((__pyx_CyFunctionObject *) (f))->defaults_getter = (g)


typedef struct {
#if CYTHON_COMPILING_IN_LIMITED_API
    PyObject_HEAD
    // We can't "inherit" from func, but we can use it as a data store
    PyObject *func;
#elif PY_VERSION_HEX < 0x030900B1
    PyCFunctionObject func;
#else
    // PEP-573: PyCFunctionObject + mm_class
    PyCMethodObject func;
#endif
#if CYTHON_COMPILING_IN_LIMITED_API && CYTHON_METH_FASTCALL
    __pyx_vectorcallfunc func_vectorcall;
#endif
#if CYTHON_COMPILING_IN_LIMITED_API
    PyObject *func_weakreflist;
#endif
#if PY_VERSION_HEX < 0x030C0000 || CYTHON_COMPILING_IN_LIMITED_API
    PyObject *func_dict;
#endif
    PyObject *func_name;
    PyObject *func_qualname;
    PyObject *func_doc;
    PyObject *func_globals;
    PyObject *func_code;
    PyObject *func_closure;
#if PY_VERSION_HEX < 0x030900B1 || CYTHON_COMPILING_IN_LIMITED_API
    // No-args super() class cell
    PyObject *func_classobj;
#endif
    // Dynamic default args and annotations
    PyObject *defaults;
    int flags;

    // Defaults info
    PyObject *defaults_tuple;   /* Const defaults tuple */
    PyObject *defaults_kwdict;  /* Const kwonly defaults dict */
    PyObject *(*defaults_getter)(PyObject *);
    PyObject *func_annotations; /* function annotations dict */

    // Coroutine marker
    PyObject *func_is_coroutine;
} __pyx_CyFunctionObject;

#undef __Pyx_CyOrPyCFunction_Check
#define __Pyx_CyFunction_Check(obj)  __Pyx_TypeCheck(obj, CGLOBAL(__pyx_CyFunctionType))
#define __Pyx_CyOrPyCFunction_Check(obj)  __Pyx_TypeCheck2(obj, CGLOBAL(__pyx_CyFunctionType), &PyCFunction_Type)
#define __Pyx_CyFunction_CheckExact(obj)  __Pyx_IS_TYPE(obj, CGLOBAL(__pyx_CyFunctionType))
static CYTHON_INLINE int __Pyx__IsSameCyOrCFunction(PyObject *func, void (*cfunc)(void));/*proto*/
#undef __Pyx_IsSameCFunction
#define __Pyx_IsSameCFunction(func, cfunc)   __Pyx__IsSameCyOrCFunction(func, cfunc)

static PyObject *__Pyx_CyFunction_Init(__pyx_CyFunctionObject* op, PyMethodDef *ml,
                                      int flags, PyObject* qualname,
                                      PyObject *closure,
                                      PyObject *module, PyObject *globals,
                                      PyObject* code);

static CYTHON_INLINE void __Pyx__CyFunction_SetClassObj(__pyx_CyFunctionObject* f, PyObject* classobj);
static CYTHON_INLINE PyObject *__Pyx_CyFunction_InitDefaults(PyObject *func,
                                                         PyTypeObject *defaults_type);
static CYTHON_INLINE void __Pyx_CyFunction_SetDefaultsTuple(PyObject *m,
                                                            PyObject *tuple);
static CYTHON_INLINE void __Pyx_CyFunction_SetDefaultsKwDict(PyObject *m,
                                                             PyObject *dict);
static CYTHON_INLINE void __Pyx_CyFunction_SetAnnotationsDict(PyObject *m,
                                                              PyObject *dict);


static int __pyx_CyFunction_init(PyObject *module);

#if CYTHON_METH_FASTCALL
static PyObject * __Pyx_CyFunction_Vectorcall_NOARGS(PyObject *func, PyObject *const *args, size_t nargsf, PyObject *kwnames);
static PyObject * __Pyx_CyFunction_Vectorcall_O(PyObject *func, PyObject *const *args, size_t nargsf, PyObject *kwnames);
static PyObject * __Pyx_CyFunction_Vectorcall_FASTCALL_KEYWORDS(PyObject *func, PyObject *const *args, size_t nargsf, PyObject *kwnames);
static PyObject * __Pyx_CyFunction_Vectorcall_FASTCALL_KEYWORDS_METHOD(PyObject *func, PyObject *const *args, size_t nargsf, PyObject *kwnames);
#if CYTHON_COMPILING_IN_LIMITED_API
#define __Pyx_CyFunction_func_vectorcall(f) (((__pyx_CyFunctionObject*)f)->func_vectorcall)
#else
#define __Pyx_CyFunction_func_vectorcall(f) (((PyCFunctionObject*)f)->vectorcall)
#endif
#endif

//////////////////// CythonFunctionShared ////////////////////
//@requires: CommonStructures.c::FetchCommonType
//@requires: CommonStructures.c::CommonTypesMetaclass
//@requires: ObjectHandling.c::PyMethodNew
//@requires: ObjectHandling.c::PyVectorcallFastCallDict
//@requires: ModuleSetupCode.c::IncludeStructmemberH
//@requires: ObjectHandling.c::PyObjectGetAttrStr
//@requires: ObjectHandling.c::CachedMethodType
//@requires: ObjectHandling.c::PyObjectCallOneArg
//@requires: ExtensionTypes.c::CallTypeTraverse
//@requires: Synchronization.c::CriticalSections
//@substitute: naming

#if CYTHON_COMPILING_IN_LIMITED_API
static CYTHON_INLINE int __Pyx__IsSameCyOrCFunctionNoMethod(PyObject *func, void (*cfunc)(void)) {
    if (__Pyx_CyFunction_Check(func)) {
        return PyCFunction_GetFunction(((__pyx_CyFunctionObject*)func)->func) == (PyCFunction) cfunc;
    } else if (PyCFunction_Check(func)) {
        return PyCFunction_GetFunction(func) == (PyCFunction) cfunc;
    }
    return 0;
}

static CYTHON_INLINE int __Pyx__IsSameCyOrCFunction(PyObject *func, void (*cfunc)(void)) {
    if ((PyObject*)Py_TYPE(func) == CGLOBAL(__Pyx_CachedMethodType)) {
        int result;
        PyObject *newFunc = PyObject_GetAttr(func, PYIDENT("__func__"));
        if (unlikely(!newFunc)) {
            PyErr_Clear(); // It's only an optimization, so don't throw an error
            return 0;
        }
        result = __Pyx__IsSameCyOrCFunctionNoMethod(newFunc, cfunc);
        Py_DECREF(newFunc);
        return result;
    }
    return __Pyx__IsSameCyOrCFunctionNoMethod(func, cfunc);
}
#else
static CYTHON_INLINE int __Pyx__IsSameCyOrCFunction(PyObject *func, void (*cfunc)(void)) {
    if (PyMethod_Check(func)) {
        func = PyMethod_GET_FUNCTION(func);
    }
    return __Pyx_CyOrPyCFunction_Check(func) && __Pyx_CyOrPyCFunction_GET_FUNCTION(func) == (PyCFunction) cfunc;
}
#endif

static CYTHON_INLINE void __Pyx__CyFunction_SetClassObj(__pyx_CyFunctionObject* f, PyObject* classobj) {
#if PY_VERSION_HEX < 0x030900B1 || CYTHON_COMPILING_IN_LIMITED_API
    __Pyx_Py_XDECREF_SET(
        __Pyx_CyFunction_GetClassObj(f),
            ((classobj) ? __Pyx_NewRef(classobj) : NULL));
#else
    __Pyx_Py_XDECREF_SET(
        // assigning to "mm_class", which is a "PyTypeObject*"
        ((PyCMethodObject *) (f))->mm_class,
        (PyTypeObject*)((classobj) ? __Pyx_NewRef(classobj) : NULL));
#endif
}

static PyObject *
__Pyx_CyFunction_get_doc_locked(__pyx_CyFunctionObject *op)
{
    if (unlikely(op->func_doc == NULL)) {
#if CYTHON_COMPILING_IN_LIMITED_API
        op->func_doc = PyObject_GetAttrString(op->func, "__doc__");
        if (unlikely(!op->func_doc)) return NULL;
#else
        if (((PyCFunctionObject*)op)->m_ml->ml_doc) {
            op->func_doc = PyUnicode_FromString(((PyCFunctionObject*)op)->m_ml->ml_doc);
            if (unlikely(op->func_doc == NULL))
                return NULL;
        } else {
            Py_INCREF(Py_None);
            return Py_None;
        }
#endif   /* CYTHON_COMPILING_IN_LIMITED_API */
    }
    Py_INCREF(op->func_doc);
    return op->func_doc;
}

static PyObject *
__Pyx_CyFunction_get_doc(__pyx_CyFunctionObject *op, void *closure) {
    PyObject *result;
    CYTHON_UNUSED_VAR(closure);
    __Pyx_BEGIN_CRITICAL_SECTION(op);
    result = __Pyx_CyFunction_get_doc_locked(op);
    __Pyx_END_CRITICAL_SECTION();
    return result;
}

static int
__Pyx_CyFunction_set_doc(__pyx_CyFunctionObject *op, PyObject *value, void *context)
{
    CYTHON_UNUSED_VAR(context);
    if (value == NULL) {
        // Mark as deleted
        value = Py_None;
    }
    Py_INCREF(value);
    __Pyx_BEGIN_CRITICAL_SECTION(op);
    __Pyx_Py_XDECREF_SET(op->func_doc, value);
    __Pyx_END_CRITICAL_SECTION();
    return 0;
}

static PyObject *
__Pyx_CyFunction_get_name_locked(__pyx_CyFunctionObject *op)
{
    if (unlikely(op->func_name == NULL)) {
#if CYTHON_COMPILING_IN_LIMITED_API
        op->func_name = PyObject_GetAttrString(op->func, "__name__");
#else
        op->func_name = PyUnicode_InternFromString(((PyCFunctionObject*)op)->m_ml->ml_name);
#endif  /* CYTHON_COMPILING_IN_LIMITED_API */
        if (unlikely(op->func_name == NULL))
            return NULL;
    }
    Py_INCREF(op->func_name);
    return op->func_name;
}

static PyObject *
__Pyx_CyFunction_get_name(__pyx_CyFunctionObject *op, void *context)
{
    PyObject *result = NULL;
    CYTHON_UNUSED_VAR(context);
    __Pyx_BEGIN_CRITICAL_SECTION(op);
    result = __Pyx_CyFunction_get_name_locked(op);
    __Pyx_END_CRITICAL_SECTION();
    return result;
}

static int
__Pyx_CyFunction_set_name(__pyx_CyFunctionObject *op, PyObject *value, void *context)
{
    CYTHON_UNUSED_VAR(context);
    if (unlikely(value == NULL || !PyUnicode_Check(value))) {
        PyErr_SetString(PyExc_TypeError,
                        "__name__ must be set to a string object");
        return -1;
    }
    Py_INCREF(value);
    __Pyx_BEGIN_CRITICAL_SECTION(op);
    __Pyx_Py_XDECREF_SET(op->func_name, value);
    __Pyx_END_CRITICAL_SECTION();
    return 0;
}

static PyObject *
__Pyx_CyFunction_get_qualname(__pyx_CyFunctionObject *op, void *context)
{
    CYTHON_UNUSED_VAR(context);
    PyObject *result;
    __Pyx_BEGIN_CRITICAL_SECTION(op);
    Py_INCREF(op->func_qualname);
    result = op->func_qualname;
    __Pyx_END_CRITICAL_SECTION();
    return result;
}

static int
__Pyx_CyFunction_set_qualname(__pyx_CyFunctionObject *op, PyObject *value, void *context)
{
    CYTHON_UNUSED_VAR(context);
    if (unlikely(value == NULL || !PyUnicode_Check(value))) {
        PyErr_SetString(PyExc_TypeError,
                        "__qualname__ must be set to a string object");
        return -1;
    }
    Py_INCREF(value);
    __Pyx_BEGIN_CRITICAL_SECTION(op);
    __Pyx_Py_XDECREF_SET(op->func_qualname, value);
    __Pyx_END_CRITICAL_SECTION();
    return 0;
}

#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000
static PyObject *
__Pyx_CyFunction_get_dict(__pyx_CyFunctionObject *op, void *context)
{
    CYTHON_UNUSED_VAR(context);
    // PyObject_GenericGetDict is not defined
    if (unlikely(op->func_dict == NULL)) {
        op->func_dict = PyDict_New();
        if (unlikely(op->func_dict == NULL))
            return NULL;
    }
    Py_INCREF(op->func_dict);
    return op->func_dict;
}
#endif

static PyObject *
__Pyx_CyFunction_get_globals(__pyx_CyFunctionObject *op, void *context)
{
    CYTHON_UNUSED_VAR(context);
    // Globals is read-only so no critical sections
    Py_INCREF(op->func_globals);
    return op->func_globals;
}

static PyObject *
__Pyx_CyFunction_get_closure(__pyx_CyFunctionObject *op, void *context)
{
    CYTHON_UNUSED_VAR(op);
    CYTHON_UNUSED_VAR(context);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
__Pyx_CyFunction_get_code(__pyx_CyFunctionObject *op, void *context)
{
    // code is read-only so no critical sections
    PyObject* result = (op->func_code) ? op->func_code : Py_None;
    CYTHON_UNUSED_VAR(context);
    Py_INCREF(result);
    return result;
}

static int
__Pyx_CyFunction_init_defaults(__pyx_CyFunctionObject *op) {
    int result = 0;
    PyObject *res = op->defaults_getter((PyObject *) op);
    if (unlikely(!res))
        return -1;

    // Cache result
    #if CYTHON_ASSUME_SAFE_MACROS && !CYTHON_AVOID_BORROWED_REFS
    op->defaults_tuple = PyTuple_GET_ITEM(res, 0);
    Py_INCREF(op->defaults_tuple);
    op->defaults_kwdict = PyTuple_GET_ITEM(res, 1);
    Py_INCREF(op->defaults_kwdict);
    #else
    op->defaults_tuple = __Pyx_PySequence_ITEM(res, 0);
    if (unlikely(!op->defaults_tuple)) result = -1;
    else {
        op->defaults_kwdict = __Pyx_PySequence_ITEM(res, 1);
        if (unlikely(!op->defaults_kwdict)) result = -1;
    }
    #endif
    Py_DECREF(res);
    return result;
}

static int
__Pyx_CyFunction_set_defaults(__pyx_CyFunctionObject *op, PyObject* value, void *context) {
    CYTHON_UNUSED_VAR(context);
    if (!value) {
        // del => explicit None to prevent rebuilding
        value = Py_None;
    } else if (unlikely(value != Py_None && !PyTuple_Check(value))) {
        PyErr_SetString(PyExc_TypeError,
                        "__defaults__ must be set to a tuple object");
        return -1;
    }
    PyErr_WarnEx(PyExc_RuntimeWarning, "changes to cyfunction.__defaults__ will not "
                 "currently affect the values used in function calls", 1);
    Py_INCREF(value);
    __Pyx_BEGIN_CRITICAL_SECTION(op);
    __Pyx_Py_XDECREF_SET(op->defaults_tuple, value);
    __Pyx_END_CRITICAL_SECTION();
    return 0;
}

static PyObject *
__Pyx_CyFunction_get_defaults_locked(__pyx_CyFunctionObject *op) {
    PyObject* result = op->defaults_tuple;
    if (unlikely(!result)) {
        if (op->defaults_getter) {
            if (unlikely(__Pyx_CyFunction_init_defaults(op) < 0)) return NULL;
            result = op->defaults_tuple;
        } else {
            result = Py_None;
        }
    }
    Py_INCREF(result);
    return result;
}

static PyObject *
__Pyx_CyFunction_get_defaults(__pyx_CyFunctionObject *op, void *context) {
    PyObject* result = NULL;
    CYTHON_UNUSED_VAR(context);
    __Pyx_BEGIN_CRITICAL_SECTION(op);
    result = __Pyx_CyFunction_get_defaults_locked(op);
    __Pyx_END_CRITICAL_SECTION();
    return result;
}

static int
__Pyx_CyFunction_set_kwdefaults(__pyx_CyFunctionObject *op, PyObject* value, void *context) {
    CYTHON_UNUSED_VAR(context);
    if (!value) {
        // del => explicit None to prevent rebuilding
        value = Py_None;
    } else if (unlikely(value != Py_None && !PyDict_Check(value))) {
        PyErr_SetString(PyExc_TypeError,
                        "__kwdefaults__ must be set to a dict object");
        return -1;
    }
    PyErr_WarnEx(PyExc_RuntimeWarning, "changes to cyfunction.__kwdefaults__ will not "
                 "currently affect the values used in function calls", 1);
    Py_INCREF(value);
    __Pyx_BEGIN_CRITICAL_SECTION(op);
    __Pyx_Py_XDECREF_SET(op->defaults_kwdict, value);
    __Pyx_END_CRITICAL_SECTION();
    return 0;
}

static PyObject *
__Pyx_CyFunction_get_kwdefaults_locked(__pyx_CyFunctionObject *op) {
    PyObject* result = op->defaults_kwdict;
    if (unlikely(!result)) {
        if (op->defaults_getter) {
            if (unlikely(__Pyx_CyFunction_init_defaults(op) < 0)) return NULL;
            result = op->defaults_kwdict;
        } else {
            result = Py_None;
        }
    }
    Py_INCREF(result);
    return result;
}

static PyObject *
__Pyx_CyFunction_get_kwdefaults(__pyx_CyFunctionObject *op, void *context) {
    PyObject* result;
    CYTHON_UNUSED_VAR(context);
    __Pyx_BEGIN_CRITICAL_SECTION(op);
    result = __Pyx_CyFunction_get_kwdefaults_locked(op);
    __Pyx_END_CRITICAL_SECTION();
    return result;
}

static int __Pyx_CyFunction_set_annotate_in_dict_if_exists(PyObject *op_in, PyObject *value); /*proto*/

static int
__Pyx_CyFunction_set_annotations(__pyx_CyFunctionObject *op, PyObject* value, void *context) {
    CYTHON_UNUSED_VAR(context);
    if (!value || value == Py_None) {
        value = NULL;
    } else if (unlikely(!PyDict_Check(value))) {
        PyErr_SetString(PyExc_TypeError,
                        "__annotations__ must be set to a dict object");
        return -1;
    }
    Py_XINCREF(value);
    __Pyx_BEGIN_CRITICAL_SECTION(op);
    __Pyx_Py_XDECREF_SET(op->func_annotations, value);
    __Pyx_END_CRITICAL_SECTION();
    if (unlikely(__Pyx_CyFunction_set_annotate_in_dict_if_exists((PyObject*) op, Py_None) < 0)) return -1;
    return 0;
}

static int
__Pyx_CyFunction_get_dict_if_exists(PyObject *op_in, PyObject **dict) {
    /* Return 1 if the function dict exists, 0 otherwise.  This cannot fail:
     * _PyObject_GetDictPtr() may clear errors internally, but never reports them. */
#if CYTHON_COMPILING_IN_PYPY
    *dict = PyObject_GenericGetDict(op_in, NULL);
#elif CYTHON_COMPILING_IN_LIMITED_API || PY_VERSION_HEX < 0x030C0000
    *dict = ((__pyx_CyFunctionObject*) op_in)->func_dict;
#else
    PyObject **dictptr = _PyObject_GetDictPtr(op_in);
    *dict = likely(dictptr) ? *dictptr : NULL;
#endif
    return *dict ? 1 : 0;
}

static int
__Pyx_CyFunction_get_annotate_from_dict_if_exists(PyObject *op_in, PyObject **annotate) {
    PyObject *dict;
    int dict_found;
    *annotate = NULL;
    dict_found = __Pyx_CyFunction_get_dict_if_exists(op_in, &dict);
    if (!dict_found) return 0;
    return __Pyx_PyDict_GetItemRef(dict, PYIDENT("__annotate__"), annotate);
}

static int
__Pyx_CyFunction_set_annotate_in_dict_if_exists(PyObject *op_in, PyObject *value) {
    PyObject *dict;
    int dict_found;
    dict_found = __Pyx_CyFunction_get_dict_if_exists(op_in, &dict);
    if (!dict_found) return 0;
    return PyDict_SetItem(dict, PYIDENT("__annotate__"), value);
}

static int
__Pyx_CyFunction_set_annotate_in_dict(PyObject *op_in, PyObject *value) {
    PyObject *dict;
    int result;
#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000
    dict = __Pyx_CyFunction_get_dict((__pyx_CyFunctionObject*) op_in, NULL);
#else
    dict = PyObject_GenericGetDict(op_in, NULL);
#endif
    if (unlikely(!dict)) return -1;
    result = PyDict_SetItem(dict, PYIDENT("__annotate__"), value);
    Py_DECREF(dict);
    return result;
}

static PyObject *
__Pyx_CyFunction_get_annotations_locked(__pyx_CyFunctionObject *op) {
    PyObject* result = op->func_annotations;
    if (unlikely(!result)) {
        result = PyDict_New();
        if (unlikely(!result)) return NULL;
        op->func_annotations = result;
    }
    Py_INCREF(result);
    return result;
}

static PyObject *
__Pyx_CyFunction_get_annotations(PyObject *op_in, void *context) {
    PyObject *annotate = NULL;
    PyObject *result = NULL;
    __pyx_CyFunctionObject *op = (__pyx_CyFunctionObject*) op_in;
    CYTHON_UNUSED_VAR(context);
    __Pyx_BEGIN_CRITICAL_SECTION(op_in);
    result = __Pyx_XNewRef(op->func_annotations);
    __Pyx_END_CRITICAL_SECTION();
    if (result) return result;

    if (unlikely(__Pyx_CyFunction_get_annotate_from_dict_if_exists(op_in, &annotate) < 0)) return NULL;
    if (!annotate || annotate == Py_None) {
        Py_XDECREF(annotate);
        __Pyx_BEGIN_CRITICAL_SECTION(op_in);
        result = __Pyx_CyFunction_get_annotations_locked(op);
        __Pyx_END_CRITICAL_SECTION();
        return result;
    }

    PyObject *format = PyLong_FromLong(1L);  // annotationlib.Format.VALUE
    if (likely(format)) {
        result = __Pyx_PyObject_CallOneArg(annotate, format);
        Py_DECREF(format);
    }
    Py_DECREF(annotate);
    if (unlikely(!result)) return NULL;
    if (unlikely(!PyDict_Check(result))) {
        PyErr_SetString(PyExc_TypeError, "__annotate__ must return a dict");
        Py_DECREF(result);
        return NULL;
    }
    __Pyx_BEGIN_CRITICAL_SECTION(op_in);
    __Pyx_Py_XDECREF_SET(op->func_annotations, __Pyx_NewRef(result));
    __Pyx_END_CRITICAL_SECTION();
    return result;
}

static PyObject *__Pyx_CyFunction_annotate_impl(PyObject *self, PyObject *args) {
    // The signature of this function is 0 or 1 args, where arg should be an int (if passed).
    // However, for this basic placeholder implementation, we don't care and don't check it.
    CYTHON_UNUSED_VAR(args);
    if (unlikely(!self)) {
        PyErr_SetString(PyExc_SystemError, "cython __annotate__ called without 'self' argument");
    }
    Py_XINCREF(self);
    return self;
}

static PyMethodDef __Pyx_CyFunction_annotate_method = {
    "__annotate__",
    (PyCFunction)(void (*)(void))__Pyx_CyFunction_annotate_impl,
    METH_VARARGS,
    "Placeholder __annotate__ function to allow 'functools.wraps' to work "
    "on Cython functions."
};

// We don't yet support PEP649 properly, so __annotate__ is
// implemented in a minimal way, just so that functools.wraps works.
static PyObject *
__Pyx_CyFunction_get_annotate(PyObject *op_in, void *context) {
    PyObject *annotate = NULL;
    CYTHON_UNUSED_VAR(context);
    if (unlikely(__Pyx_CyFunction_get_annotate_from_dict_if_exists(op_in, &annotate) < 0)) return NULL;
    if (annotate) return annotate;

    PyObject *annotations = __Pyx_CyFunction_get_annotations(op_in, NULL);
    if (unlikely(!annotations)) return NULL;
    PyObject *method = PyCFunction_New(
        &__Pyx_CyFunction_annotate_method,
        annotations);
    Py_DECREF(annotations);
    return method;
}

static int
__Pyx_CyFunction_set_annotate(PyObject *op_in, PyObject* value, void *context) {
    CYTHON_UNUSED_VAR(context);
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "__annotate__ cannot be deleted");
        return -1;
    }
    if (unlikely(value != Py_None && !PyCallable_Check(value))) {
        PyErr_SetString(PyExc_TypeError, "__annotate__ must be callable or None");
        return -1;
    }
    if (value != Py_None) {
        __pyx_CyFunctionObject *op = (__pyx_CyFunctionObject*) op_in;
        __Pyx_BEGIN_CRITICAL_SECTION(op_in);
        Py_CLEAR(op->func_annotations);
        __Pyx_END_CRITICAL_SECTION();
    }
    return __Pyx_CyFunction_set_annotate_in_dict(op_in, value);
}

static PyObject *
__Pyx_CyFunction_get_is_coroutine_value(__pyx_CyFunctionObject *op) {
    int is_coroutine = op->flags & __Pyx_CYFUNCTION_COROUTINE;
    if (is_coroutine) {
        PyObject *is_coroutine_value, *module, *fromlist, *marker = PYIDENT("_is_coroutine");
        fromlist = PyList_New(1);
        if (unlikely(!fromlist)) return NULL;
        Py_INCREF(marker);
#if CYTHON_ASSUME_SAFE_MACROS
        PyList_SET_ITEM(fromlist, 0, marker);
#else
        if (unlikely(PyList_SetItem(fromlist, 0, marker) < 0)) {
            Py_DECREF(fromlist);
            return NULL;
        }
#endif
        module = PyImport_ImportModuleLevelObject(PYIDENT("asyncio.coroutines"), NULL, NULL, fromlist, 0);
        Py_DECREF(fromlist);
        if (unlikely(!module)) goto ignore;
        is_coroutine_value = __Pyx_PyObject_GetAttrStr(module, marker);
        Py_DECREF(module);
        if (likely(is_coroutine_value)) {
            return is_coroutine_value;
        }
ignore:
        PyErr_Clear();
    }

    return __Pyx_PyBool_FromLong(is_coroutine);
}

static PyObject *
__Pyx_CyFunction_get_is_coroutine(__pyx_CyFunctionObject *op, void *context) {
    PyObject *result;
    CYTHON_UNUSED_VAR(context);
    if (op->func_is_coroutine) {
        return __Pyx_NewRef(op->func_is_coroutine);
    }

    result = __Pyx_CyFunction_get_is_coroutine_value(op);
    if (unlikely(!result))
        return NULL;

    __Pyx_BEGIN_CRITICAL_SECTION(op);
    // Guard against concurrent initialisation.
    if (op->func_is_coroutine) {
        Py_DECREF(result);
        result = __Pyx_NewRef(op->func_is_coroutine);
    } else {
        op->func_is_coroutine = __Pyx_NewRef(result);
    }
    __Pyx_END_CRITICAL_SECTION();
    return result;
}

//static PyObject *
//__Pyx_CyFunction_get_signature(__pyx_CyFunctionObject *op, void *context) {
//    PyObject *inspect_module, *signature_class, *signature;
//    CYTHON_UNUSED_VAR(context);
//    // from inspect import Signature
//    inspect_module = PyImport_ImportModuleLevelObject(PYIDENT("inspect"), NULL, NULL, NULL, 0);
//    if (unlikely(!inspect_module))
//        goto bad;
//    signature_class = __Pyx_PyObject_GetAttrStr(inspect_module, PYIDENT("Signature"));
//    Py_DECREF(inspect_module);
//    if (unlikely(!signature_class))
//        goto bad;
//    // return Signature.from_function(op)
//    signature = PyObject_CallMethodObjArgs(signature_class, PYIDENT("from_function"), op, NULL);
//    Py_DECREF(signature_class);
//    if (likely(signature))
//        return signature;
//bad:
//    // make sure we raise an AttributeError from this property on any errors
//    if (!PyErr_ExceptionMatches(PyExc_AttributeError))
//        PyErr_SetString(PyExc_AttributeError, "failed to calculate __signature__");
//    return NULL;
//}

static void __Pyx_CyFunction_raise_argument_count_error(__pyx_CyFunctionObject *func, const char* message, Py_ssize_t size) {
#if CYTHON_COMPILING_IN_LIMITED_API
    PyObject *py_name = __Pyx_CyFunction_get_name(func, NULL);
    if (!py_name) return;
    PyErr_Format(PyExc_TypeError,
        "%.200S() %s (%" CYTHON_FORMAT_SSIZE_T "d given)",
        py_name, message, size);
    Py_DECREF(py_name);
#else
    const char* name = ((PyCFunctionObject*)func)->m_ml->ml_name;
    PyErr_Format(PyExc_TypeError,
        "%.200s() %s (%" CYTHON_FORMAT_SSIZE_T "d given)",
        name, message, size);
#endif
}

static void __Pyx_CyFunction_raise_type_error(__pyx_CyFunctionObject *func, const char* message) {
#if CYTHON_COMPILING_IN_LIMITED_API
    PyObject *py_name = __Pyx_CyFunction_get_name(func, NULL);
    if (!py_name) return;
    PyErr_Format(PyExc_TypeError,
        "%.200S() %s",
        py_name, message);
    Py_DECREF(py_name);
#else
    const char* name = ((PyCFunctionObject*)func)->m_ml->ml_name;
    PyErr_Format(PyExc_TypeError,
        "%.200s() %s",
        name, message);
#endif
}


#if CYTHON_COMPILING_IN_LIMITED_API
static PyObject *
__Pyx_CyFunction_get_module(__pyx_CyFunctionObject *op, void *context) {
    CYTHON_UNUSED_VAR(context);
    return PyObject_GetAttrString(op->func, "__module__");
}

static int
__Pyx_CyFunction_set_module(__pyx_CyFunctionObject *op, PyObject* value, void *context) {
    CYTHON_UNUSED_VAR(context);
    return PyObject_SetAttrString(op->func, "__module__", value);
}
#endif

static PyGetSetDef __pyx_CyFunction_getsets[] = {
    {"func_doc", (getter)__Pyx_CyFunction_get_doc, (setter)__Pyx_CyFunction_set_doc, 0, 0},
    {"__doc__",  (getter)__Pyx_CyFunction_get_doc, (setter)__Pyx_CyFunction_set_doc, 0, 0},
    {"func_name", (getter)__Pyx_CyFunction_get_name, (setter)__Pyx_CyFunction_set_name, 0, 0},
    {"__name__", (getter)__Pyx_CyFunction_get_name, (setter)__Pyx_CyFunction_set_name, 0, 0},
    {"__qualname__", (getter)__Pyx_CyFunction_get_qualname, (setter)__Pyx_CyFunction_set_qualname, 0, 0},
#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000
    {"func_dict", (getter)__Pyx_CyFunction_get_dict, (setter)PyObject_GenericSetDict, 0, 0},
    {"__dict__", (getter)__Pyx_CyFunction_get_dict, (setter)PyObject_GenericSetDict, 0, 0},
#else
    {"func_dict", (getter)PyObject_GenericGetDict, (setter)PyObject_GenericSetDict, 0, 0},
    {"__dict__", (getter)PyObject_GenericGetDict, (setter)PyObject_GenericSetDict, 0, 0},
#endif
    {"func_globals", (getter)__Pyx_CyFunction_get_globals, 0, 0, 0},
    {"__globals__", (getter)__Pyx_CyFunction_get_globals, 0, 0, 0},
    {"func_closure", (getter)__Pyx_CyFunction_get_closure, 0, 0, 0},
    {"__closure__", (getter)__Pyx_CyFunction_get_closure, 0, 0, 0},
    {"func_code", (getter)__Pyx_CyFunction_get_code, 0, 0, 0},
    {"__code__", (getter)__Pyx_CyFunction_get_code, 0, 0, 0},
    {"func_defaults", (getter)__Pyx_CyFunction_get_defaults, (setter)__Pyx_CyFunction_set_defaults, 0, 0},
    {"__defaults__", (getter)__Pyx_CyFunction_get_defaults, (setter)__Pyx_CyFunction_set_defaults, 0, 0},
    {"__kwdefaults__", (getter)__Pyx_CyFunction_get_kwdefaults, (setter)__Pyx_CyFunction_set_kwdefaults, 0, 0},
    {"__annotations__", (getter)__Pyx_CyFunction_get_annotations, (setter)__Pyx_CyFunction_set_annotations, 0, 0},
    {"__annotate__", (getter)__Pyx_CyFunction_get_annotate, (setter)__Pyx_CyFunction_set_annotate, 0, 0},
    {"_is_coroutine", (getter)__Pyx_CyFunction_get_is_coroutine, 0, 0, 0},
//    {"__signature__", (getter)__Pyx_CyFunction_get_signature, 0, 0, 0},
#if CYTHON_COMPILING_IN_LIMITED_API
    {"__module__", (getter)__Pyx_CyFunction_get_module, (setter)__Pyx_CyFunction_set_module, 0, 0},
#endif
    {0, 0, 0, 0, 0}
};

static PyMemberDef __pyx_CyFunction_members[] = {
#if !CYTHON_COMPILING_IN_LIMITED_API
    {"__module__", T_OBJECT, offsetof(PyCFunctionObject, m_module), 0, 0},
#endif
#if PY_VERSION_HEX < 0x030C0000 || CYTHON_COMPILING_IN_LIMITED_API
    {"__dictoffset__", T_PYSSIZET, offsetof(__pyx_CyFunctionObject, func_dict), READONLY, 0},
#endif
#if CYTHON_METH_FASTCALL
#if CYTHON_COMPILING_IN_LIMITED_API
    {"__vectorcalloffset__", T_PYSSIZET, offsetof(__pyx_CyFunctionObject, func_vectorcall), READONLY, 0},
#else
    {"__vectorcalloffset__", T_PYSSIZET, offsetof(PyCFunctionObject, vectorcall), READONLY, 0},
#endif
#if CYTHON_COMPILING_IN_LIMITED_API
    {"__weaklistoffset__", T_PYSSIZET, offsetof(__pyx_CyFunctionObject, func_weakreflist), READONLY, 0},
#else
    {"__weaklistoffset__", T_PYSSIZET, offsetof(PyCFunctionObject, m_weakreflist), READONLY, 0},
#endif
#endif
    {0, 0, 0,  0, 0}
};

static PyObject *
__Pyx_CyFunction_reduce(__pyx_CyFunctionObject *m, PyObject *args)
{
    PyObject *result = NULL;
    CYTHON_UNUSED_VAR(args);
    __Pyx_BEGIN_CRITICAL_SECTION(m);
    Py_INCREF(m->func_qualname);
    result = m->func_qualname;
    __Pyx_END_CRITICAL_SECTION();
    return result;
}

static PyMethodDef __pyx_CyFunction_methods[] = {
    {"__reduce__", (PyCFunction)__Pyx_CyFunction_reduce, METH_VARARGS, 0},
    {0, 0, 0, 0}
};


#if CYTHON_COMPILING_IN_LIMITED_API
#define __Pyx_CyFunction_weakreflist(cyfunc) ((cyfunc)->func_weakreflist)
#else
#define __Pyx_CyFunction_weakreflist(cyfunc) (((PyCFunctionObject*)cyfunc)->m_weakreflist)
#endif

static PyObject *__Pyx_CyFunction_Init(__pyx_CyFunctionObject *op, PyMethodDef *ml, int flags, PyObject* qualname,
                                       PyObject *closure, PyObject *module, PyObject* globals, PyObject* code) {
#if !CYTHON_COMPILING_IN_LIMITED_API
    PyCFunctionObject *cf = (PyCFunctionObject*) op;
#endif
    if (unlikely(op == NULL))
        return NULL;
#if CYTHON_COMPILING_IN_LIMITED_API
    // Note that we end up with a circular reference to op. This isn't
    // a disaster, but in an ideal world it'd be nice to avoid it.
    op->func = PyCFunction_NewEx(ml, (PyObject*)op, module);
    if (unlikely(!op->func)) return NULL;
#endif
    op->flags = flags;
    __Pyx_CyFunction_weakreflist(op) = NULL;
#if !CYTHON_COMPILING_IN_LIMITED_API
    cf->m_ml = ml;
    cf->m_self = (PyObject *) op;
#endif
    Py_XINCREF(closure);
    op->func_closure = closure;
#if !CYTHON_COMPILING_IN_LIMITED_API
    Py_XINCREF(module);
    cf->m_module = module;
#endif
#if PY_VERSION_HEX < 0x030C0000 || CYTHON_COMPILING_IN_LIMITED_API
    op->func_dict = NULL;
#endif
    op->func_name = NULL;
    Py_INCREF(qualname);
    op->func_qualname = qualname;
    op->func_doc = NULL;
#if PY_VERSION_HEX < 0x030900B1 || CYTHON_COMPILING_IN_LIMITED_API
    op->func_classobj = NULL;
#else
    ((PyCMethodObject*)op)->mm_class = NULL;
#endif
    op->func_globals = globals;
    Py_INCREF(op->func_globals);
    Py_XINCREF(code);
    op->func_code = code;
    // Dynamic Default args
    op->defaults = NULL;
    op->defaults_tuple = NULL;
    op->defaults_kwdict = NULL;
    op->defaults_getter = NULL;
    op->func_annotations = NULL;
    op->func_is_coroutine = NULL;
#if CYTHON_METH_FASTCALL
    switch (ml->ml_flags & (METH_VARARGS | METH_FASTCALL | METH_NOARGS | METH_O | METH_KEYWORDS | METH_METHOD)) {
    case METH_NOARGS:
        __Pyx_CyFunction_func_vectorcall(op) = __Pyx_CyFunction_Vectorcall_NOARGS;
        break;
    case METH_O:
        __Pyx_CyFunction_func_vectorcall(op) = __Pyx_CyFunction_Vectorcall_O;
        break;
    // case METH_FASTCALL is not used
    case METH_METHOD | METH_FASTCALL | METH_KEYWORDS:
        __Pyx_CyFunction_func_vectorcall(op) = __Pyx_CyFunction_Vectorcall_FASTCALL_KEYWORDS_METHOD;
        break;
    case METH_FASTCALL | METH_KEYWORDS:
        __Pyx_CyFunction_func_vectorcall(op) = __Pyx_CyFunction_Vectorcall_FASTCALL_KEYWORDS;
        break;
    // case METH_VARARGS is not used
    case METH_VARARGS | METH_KEYWORDS:
        __Pyx_CyFunction_func_vectorcall(op) = NULL;
        break;
    default:
        PyErr_SetString(PyExc_SystemError, "Bad call flags for CyFunction");
        Py_DECREF(op);
        return NULL;
    }
#endif
    return (PyObject *) op;
}

static int
__Pyx_CyFunction_clear(__pyx_CyFunctionObject *m)
{
    Py_CLEAR(m->func_closure);
#if CYTHON_COMPILING_IN_LIMITED_API
    Py_CLEAR(m->func);
#else
    Py_CLEAR(((PyCFunctionObject*)m)->m_module);
#endif
#if PY_VERSION_HEX < 0x030C0000 || CYTHON_COMPILING_IN_LIMITED_API
    Py_CLEAR(m->func_dict);
#elif PY_VERSION_HEX < 0x030d0000
    _PyObject_ClearManagedDict((PyObject*)m);
#else
    PyObject_ClearManagedDict((PyObject*)m);
#endif
    Py_CLEAR(m->func_name);
    Py_CLEAR(m->func_qualname);
    Py_CLEAR(m->func_doc);
    Py_CLEAR(m->func_globals);
    Py_CLEAR(m->func_code);
#if PY_VERSION_HEX < 0x030900B1 || CYTHON_COMPILING_IN_LIMITED_API
    Py_CLEAR(__Pyx_CyFunction_GetClassObj(m));
#else
    {
        PyObject *cls = (PyObject*) ((PyCMethodObject *) (m))->mm_class;
        ((PyCMethodObject *) (m))->mm_class = NULL;
        Py_XDECREF(cls);
    }
#endif
    Py_CLEAR(m->defaults_tuple);
    Py_CLEAR(m->defaults_kwdict);
    Py_CLEAR(m->func_annotations);
    Py_CLEAR(m->func_is_coroutine);

    Py_CLEAR(m->defaults);

    return 0;
}

static void __Pyx__CyFunction_dealloc(__pyx_CyFunctionObject *m)
{
    if (__Pyx_CyFunction_weakreflist(m) != NULL)
        PyObject_ClearWeakRefs((PyObject *) m);
    __Pyx_CyFunction_clear(m);
    __Pyx_PyHeapTypeObject_GC_Del(m);
}

static void __Pyx_CyFunction_dealloc(__pyx_CyFunctionObject *m)
{
    PyObject_GC_UnTrack(m);
    __Pyx__CyFunction_dealloc(m);
}

static int __Pyx_CyFunction_traverse(__pyx_CyFunctionObject *m, visitproc visit, void *arg)
{
    {
        int e = __Pyx_call_type_traverse((PyObject*)m, 1, visit, arg);
        if (e) return e;
    }
    Py_VISIT(m->func_closure);
#if CYTHON_COMPILING_IN_LIMITED_API
    Py_VISIT(m->func);
#else
    Py_VISIT(((PyCFunctionObject*)m)->m_module);
#endif
#if PY_VERSION_HEX < 0x030C0000 || CYTHON_COMPILING_IN_LIMITED_API
    Py_VISIT(m->func_dict);
#else
    {
        int e =
#if PY_VERSION_HEX < 0x030d0000
            _PyObject_VisitManagedDict
#else
            PyObject_VisitManagedDict
#endif
                ((PyObject*)m, visit, arg);
        if (e != 0) return e;
    }
#endif
    __Pyx_VISIT_CONST(m->func_name);
    __Pyx_VISIT_CONST(m->func_qualname);
    Py_VISIT(m->func_doc);
    Py_VISIT(m->func_globals);
    // The code objects that we generate only contain plain constants and can never participate in reference cycles.
    __Pyx_VISIT_CONST(m->func_code);
    Py_VISIT(__Pyx_CyFunction_GetClassObj(m));
    Py_VISIT(m->defaults_tuple);
    Py_VISIT(m->defaults_kwdict);
    Py_VISIT(m->func_annotations);
    Py_VISIT(m->func_is_coroutine);
    Py_VISIT(m->defaults);

    return 0;
}

static PyObject*
__Pyx_CyFunction_repr(__pyx_CyFunctionObject *op)
{
    PyObject *repr;
    __Pyx_BEGIN_CRITICAL_SECTION(op);
    repr = PyUnicode_FromFormat("<cyfunction %U at %p>",
                                op->func_qualname, (void *)op);
    __Pyx_END_CRITICAL_SECTION();
    return repr;
}

static PyObject * __Pyx_CyFunction_CallMethod(PyObject *func, PyObject *self, PyObject *arg, PyObject *kw) {
    // originally copied from PyCFunction_Call() in CPython's Objects/methodobject.c
#if CYTHON_COMPILING_IN_LIMITED_API
    PyObject *f = ((__pyx_CyFunctionObject*)func)->func;
    PyCFunction meth;
    int flags;
    meth = PyCFunction_GetFunction(f);
    if (unlikely(!meth)) return NULL;
    flags = PyCFunction_GetFlags(f);
    if (unlikely(flags < 0)) return NULL;
#else
    PyCFunctionObject* f = (PyCFunctionObject*)func;
    PyCFunction meth = f->m_ml->ml_meth;
    int flags = f->m_ml->ml_flags;
#endif
    Py_ssize_t size;

    switch (flags & (METH_VARARGS | METH_KEYWORDS | METH_NOARGS | METH_O)) {
    case METH_VARARGS:
        if (likely(kw == NULL || PyDict_Size(kw) == 0))
            return (*meth)(self, arg);
        break;
    case METH_VARARGS | METH_KEYWORDS:
        return (*(PyCFunctionWithKeywords)(void(*)(void))meth)(self, arg, kw);
    case METH_NOARGS:
        if (likely(kw == NULL || PyDict_Size(kw) == 0)) {
#if CYTHON_ASSUME_SAFE_SIZE
            size = PyTuple_GET_SIZE(arg);
#else
            size = PyTuple_Size(arg);
            if (unlikely(size < 0)) return NULL;
#endif
            if (likely(size == 0))
                return (*meth)(self, NULL);
            __Pyx_CyFunction_raise_argument_count_error(
                (__pyx_CyFunctionObject*)func,
                "takes no arguments", size);
            return NULL;
        }
        break;
    case METH_O:
        if (likely(kw == NULL || PyDict_Size(kw) == 0)) {
#if CYTHON_ASSUME_SAFE_SIZE
            size = PyTuple_GET_SIZE(arg);
#else
            size = PyTuple_Size(arg);
            if (unlikely(size < 0)) return NULL;
#endif
            if (likely(size == 1)) {
                PyObject *result, *arg0;
                #if CYTHON_ASSUME_SAFE_MACROS && !CYTHON_AVOID_BORROWED_REFS
                arg0 = PyTuple_GET_ITEM(arg, 0);
                #else
                arg0 = __Pyx_PySequence_ITEM(arg, 0); if (unlikely(!arg0)) return NULL;
                #endif
                result = (*meth)(self, arg0);
                #if !(CYTHON_ASSUME_SAFE_MACROS && !CYTHON_AVOID_BORROWED_REFS)
                Py_DECREF(arg0);
                #endif
                return result;
            }
            __Pyx_CyFunction_raise_argument_count_error(
                (__pyx_CyFunctionObject*)func,
                "takes exactly one argument", size);
            return NULL;
        }
        break;
    default:
        PyErr_SetString(PyExc_SystemError, "Bad call flags for CyFunction");
        return NULL;
    }
    __Pyx_CyFunction_raise_type_error(
        (__pyx_CyFunctionObject*)func, "takes no keyword arguments");
    return NULL;
}

static CYTHON_INLINE PyObject *__Pyx_CyFunction_Call(PyObject *func, PyObject *arg, PyObject *kw) {
    PyObject *self, *result;
#if CYTHON_COMPILING_IN_LIMITED_API
    // PyCFunction_GetSelf returns a borrowed reference
    self = PyCFunction_GetSelf(((__pyx_CyFunctionObject*)func)->func);
    if (unlikely(!self) && PyErr_Occurred()) return NULL;
#else
    self = ((PyCFunctionObject*)func)->m_self;
#endif
    result = __Pyx_CyFunction_CallMethod(func, self, arg, kw);
    return result;
}

static PyObject *__Pyx_CyFunction_CallAsMethod(PyObject *func, PyObject *args, PyObject *kw) {
    PyObject *result;
    __pyx_CyFunctionObject *cyfunc = (__pyx_CyFunctionObject *) func;

#if CYTHON_METH_FASTCALL && CYTHON_VECTORCALL
    // Prefer vectorcall if available. This is not the typical case, as
    // CPython would normally use vectorcall directly instead of tp_call.
     __pyx_vectorcallfunc vc = __Pyx_CyFunction_func_vectorcall(cyfunc);
    if (vc) {
#if CYTHON_ASSUME_SAFE_MACROS && CYTHON_ASSUME_SAFE_SIZE
        return __Pyx_PyVectorcall_FastCallDict(func, vc, &PyTuple_GET_ITEM(args, 0), (size_t)PyTuple_GET_SIZE(args), kw);
#else
        // avoid unused function warning
        (void) &__Pyx_PyVectorcall_FastCallDict;
        return PyVectorcall_Call(func, args, kw);
#endif
    }
#endif

    if ((cyfunc->flags & __Pyx_CYFUNCTION_CCLASS) && !(cyfunc->flags & __Pyx_CYFUNCTION_STATICMETHOD)) {
        Py_ssize_t argc;
        PyObject *new_args;
        PyObject *self;

#if CYTHON_ASSUME_SAFE_SIZE
        argc = PyTuple_GET_SIZE(args);
#else
        argc = PyTuple_Size(args);
        if (unlikely(argc < 0)) return NULL;
#endif
        new_args = PyTuple_GetSlice(args, 1, argc);

        if (unlikely(!new_args))
            return NULL;

        self = PyTuple_GetItem(args, 0);
        if (unlikely(!self)) {
            Py_DECREF(new_args);
            PyErr_Format(PyExc_TypeError,
                         "unbound method %.200S() needs an argument",
                         cyfunc->func_qualname);
            return NULL;
        }

        result = __Pyx_CyFunction_CallMethod(func, self, new_args, kw);
        Py_DECREF(new_args);
    } else {
        result = __Pyx_CyFunction_Call(func, args, kw);
    }
    return result;
}

#if CYTHON_METH_FASTCALL && CYTHON_VECTORCALL
// Check that kwnames is empty (if you want to allow keyword arguments,
// simply pass kwnames=NULL) and figure out what to do with "self".
// Return value:
//  1: self = args[0]
//  0: self = cyfunc->func.m_self
// -1: error
static CYTHON_INLINE int __Pyx_CyFunction_Vectorcall_CheckArgs(__pyx_CyFunctionObject *cyfunc, Py_ssize_t nargs, PyObject *kwnames)
{
    int ret = 0;
    if ((cyfunc->flags & __Pyx_CYFUNCTION_CCLASS) && !(cyfunc->flags & __Pyx_CYFUNCTION_STATICMETHOD)) {
        if (unlikely(nargs < 1)) {
            __Pyx_CyFunction_raise_type_error(
                cyfunc, "needs an argument");
            return -1;
        }
        ret = 1;
    }
    if (unlikely(kwnames) && unlikely(__Pyx_PyTuple_GET_SIZE(kwnames))) {
        __Pyx_CyFunction_raise_type_error(
            cyfunc, "takes no keyword arguments");
        return -1;
    }
    return ret;
}

static PyObject * __Pyx_CyFunction_Vectorcall_NOARGS(PyObject *func, PyObject *const *args, size_t nargsf, PyObject *kwnames)
{
    __pyx_CyFunctionObject *cyfunc = (__pyx_CyFunctionObject *)func;
    Py_ssize_t nargs = PyVectorcall_NARGS(nargsf);
    PyObject *self;
#if CYTHON_COMPILING_IN_LIMITED_API
    PyCFunction meth = PyCFunction_GetFunction(cyfunc->func);
    if (unlikely(!meth)) return NULL;
#else
    PyCFunction meth = ((PyCFunctionObject*)cyfunc)->m_ml->ml_meth;
#endif

    switch (__Pyx_CyFunction_Vectorcall_CheckArgs(cyfunc, nargs, kwnames)) {
    case 1:
        self = args[0];
        args += 1;
        nargs -= 1;
        break;
    case 0:
#if CYTHON_COMPILING_IN_LIMITED_API
        // PyCFunction_GetSelf returns a borrowed reference
        self = PyCFunction_GetSelf(((__pyx_CyFunctionObject*)cyfunc)->func);
        if (unlikely(!self) && PyErr_Occurred()) return NULL;
#else
        self = ((PyCFunctionObject*)cyfunc)->m_self;
#endif
        break;
    default:
        return NULL;
    }

    if (unlikely(nargs != 0)) {
        __Pyx_CyFunction_raise_argument_count_error(
            cyfunc, "takes no arguments", nargs);
        return NULL;
    }
    return meth(self, NULL);
}

static PyObject * __Pyx_CyFunction_Vectorcall_O(PyObject *func, PyObject *const *args, size_t nargsf, PyObject *kwnames)
{
    __pyx_CyFunctionObject *cyfunc = (__pyx_CyFunctionObject *)func;
    Py_ssize_t nargs = PyVectorcall_NARGS(nargsf);
    PyObject *self;
#if CYTHON_COMPILING_IN_LIMITED_API
    PyCFunction meth = PyCFunction_GetFunction(cyfunc->func);
    if (unlikely(!meth)) return NULL;
#else
    PyCFunction meth = ((PyCFunctionObject*)cyfunc)->m_ml->ml_meth;
#endif

    switch (__Pyx_CyFunction_Vectorcall_CheckArgs(cyfunc, nargs, kwnames)) {
    case 1:
        self = args[0];
        args += 1;
        nargs -= 1;
        break;
    case 0:
#if CYTHON_COMPILING_IN_LIMITED_API
        // PyCFunction_GetSelf returns a borrowed reference
        self = PyCFunction_GetSelf(((__pyx_CyFunctionObject*)cyfunc)->func);
        if (unlikely(!self) && PyErr_Occurred()) return NULL;
#else
        self = ((PyCFunctionObject*)cyfunc)->m_self;
#endif
        break;
    default:
        return NULL;
    }

    if (unlikely(nargs != 1)) {
        __Pyx_CyFunction_raise_argument_count_error(
            cyfunc, "takes exactly one argument", nargs);
        return NULL;
    }
    return meth(self, args[0]);
}

static PyObject * __Pyx_CyFunction_Vectorcall_FASTCALL_KEYWORDS(PyObject *func, PyObject *const *args, size_t nargsf, PyObject *kwnames)
{
    __pyx_CyFunctionObject *cyfunc = (__pyx_CyFunctionObject *)func;
    Py_ssize_t nargs = PyVectorcall_NARGS(nargsf);
    PyObject *self;
#if CYTHON_COMPILING_IN_LIMITED_API
    PyCFunction meth = PyCFunction_GetFunction(cyfunc->func);
    if (unlikely(!meth)) return NULL;
#else
    PyCFunction meth = ((PyCFunctionObject*)cyfunc)->m_ml->ml_meth;
#endif

    switch (__Pyx_CyFunction_Vectorcall_CheckArgs(cyfunc, nargs, NULL)) {
    case 1:
        self = args[0];
        args += 1;
        nargs -= 1;
        break;
    case 0:
#if CYTHON_COMPILING_IN_LIMITED_API
        // PyCFunction_GetSelf returns a borrowed reference
        self = PyCFunction_GetSelf(((__pyx_CyFunctionObject*)cyfunc)->func);
        if (unlikely(!self) && PyErr_Occurred()) return NULL;
#else
        self = ((PyCFunctionObject*)cyfunc)->m_self;
#endif
        break;
    default:
        return NULL;
    }

    return ((__Pyx_PyCFunctionFastWithKeywords)(void(*)(void))meth)(self, args, nargs, kwnames);
}

static PyObject * __Pyx_CyFunction_Vectorcall_FASTCALL_KEYWORDS_METHOD(PyObject *func, PyObject *const *args, size_t nargsf, PyObject *kwnames)
{
    __pyx_CyFunctionObject *cyfunc = (__pyx_CyFunctionObject *)func;
    PyTypeObject *cls = (PyTypeObject *) __Pyx_CyFunction_GetClassObj(cyfunc);
    Py_ssize_t nargs = PyVectorcall_NARGS(nargsf);
    PyObject *self;
#if CYTHON_COMPILING_IN_LIMITED_API
    PyCFunction meth = PyCFunction_GetFunction(cyfunc->func);
    if (unlikely(!meth)) return NULL;
#else
    PyCFunction meth = ((PyCFunctionObject*)cyfunc)->m_ml->ml_meth;
#endif
    switch (__Pyx_CyFunction_Vectorcall_CheckArgs(cyfunc, nargs, NULL)) {
    case 1:
        self = args[0];
        args += 1;
        nargs -= 1;
        break;
    case 0:
#if CYTHON_COMPILING_IN_LIMITED_API
        // PyCFunction_GetSelf returns a borrowed reference
        self = PyCFunction_GetSelf(((__pyx_CyFunctionObject*)cyfunc)->func);
        if (unlikely(!self) && PyErr_Occurred()) return NULL;
#else
        self = ((PyCFunctionObject*)cyfunc)->m_self;
#endif
        break;
    default:
        return NULL;
    }

    #if PY_VERSION_HEX < 0x030e00A6
    // See https://github.com/python/cpython/pull/131135
    size_t nargs_value = (size_t) nargs;
    #else
    Py_ssize_t nargs_value = nargs;
    #endif

    return ((__Pyx_PyCMethod)(void(*)(void))meth)(self, cls, args, nargs_value, kwnames);
}
#endif

static PyType_Slot __pyx_CyFunctionType_slots[] = {
    {Py_tp_dealloc, (void *)__Pyx_CyFunction_dealloc},
    {Py_tp_repr, (void *)__Pyx_CyFunction_repr},
    {Py_tp_call, (void *)__Pyx_CyFunction_CallAsMethod},
    {Py_tp_traverse, (void *)__Pyx_CyFunction_traverse},
    {Py_tp_clear, (void *)__Pyx_CyFunction_clear},
    {Py_tp_methods, (void *)__pyx_CyFunction_methods},
    {Py_tp_members, (void *)__pyx_CyFunction_members},
    {Py_tp_getset, (void *)__pyx_CyFunction_getsets},
    {Py_tp_descr_get, (void *)__Pyx_PyMethod_New},
    {0, 0},
};

static PyType_Spec __pyx_CyFunctionType_spec = {
    __PYX_TYPE_MODULE_PREFIX "cython_function_or_method",
    sizeof(__pyx_CyFunctionObject),
    0,
#ifdef Py_TPFLAGS_METHOD_DESCRIPTOR
    Py_TPFLAGS_METHOD_DESCRIPTOR |
#endif
#if CYTHON_METH_FASTCALL
#if defined(Py_TPFLAGS_HAVE_VECTORCALL)
    Py_TPFLAGS_HAVE_VECTORCALL |
#elif defined(_Py_TPFLAGS_HAVE_VECTORCALL)
    _Py_TPFLAGS_HAVE_VECTORCALL |
#endif
#endif // CYTHON_METH_FASTCALL
#if PY_VERSION_HEX >= 0x030C0000 && !CYTHON_COMPILING_IN_LIMITED_API
    Py_TPFLAGS_MANAGED_DICT |
#endif
    Py_TPFLAGS_IMMUTABLETYPE | Py_TPFLAGS_DISALLOW_INSTANTIATION |
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    __pyx_CyFunctionType_slots
};

static int __pyx_CyFunction_init(PyObject *module) {
    $modulestatetype_cname *mstate = __Pyx_PyModule_GetState(module);
    mstate->__pyx_CyFunctionType = __Pyx_FetchCommonTypeFromSpec(
        mstate->__pyx_CommonTypesMetaclassType, module, &__pyx_CyFunctionType_spec, NULL);
    if (unlikely(mstate->__pyx_CyFunctionType == NULL)) {
        return -1;
    }
    return 0;
}

static CYTHON_INLINE PyObject *__Pyx_CyFunction_InitDefaults(PyObject *func, PyTypeObject *defaults_type) {
    __pyx_CyFunctionObject *m = (__pyx_CyFunctionObject *) func;

    m->defaults = PyObject_CallObject((PyObject*)defaults_type, NULL); // _PyObject_New(defaults_type);
    if (unlikely(!m->defaults))
        return NULL;
    return m->defaults;
}

static CYTHON_INLINE void __Pyx_CyFunction_SetDefaultsTuple(PyObject *func, PyObject *tuple) {
    __pyx_CyFunctionObject *m = (__pyx_CyFunctionObject *) func;
    m->defaults_tuple = tuple;
    Py_INCREF(tuple);
}

static CYTHON_INLINE void __Pyx_CyFunction_SetDefaultsKwDict(PyObject *func, PyObject *dict) {
    __pyx_CyFunctionObject *m = (__pyx_CyFunctionObject *) func;
    m->defaults_kwdict = dict;
    Py_INCREF(dict);
}

static CYTHON_INLINE void __Pyx_CyFunction_SetAnnotationsDict(PyObject *func, PyObject *dict) {
    __pyx_CyFunctionObject *m = (__pyx_CyFunctionObject *) func;
    m->func_annotations = dict;
    Py_INCREF(dict);
}


//////////////////// CythonFunction.proto ////////////////////

static PyObject *__Pyx_CyFunction_New(PyMethodDef *ml,
                                      int flags, PyObject* qualname,
                                      PyObject *closure,
                                      PyObject *module, PyObject *globals,
                                      PyObject* code);

//////////////////// CythonFunction ////////////////////
//@requires: CythonFunctionShared

static PyObject *__Pyx_CyFunction_New(PyMethodDef *ml, int flags, PyObject* qualname,
                                      PyObject *closure, PyObject *module, PyObject* globals, PyObject* code) {
    PyObject *op = __Pyx_CyFunction_Init(
        PyObject_GC_New(__pyx_CyFunctionObject, CGLOBAL(__pyx_CyFunctionType)),
        ml, flags, qualname, closure, module, globals, code
    );
    if (likely(op)) {
        PyObject_GC_Track(op);
    }
    return op;
}

//////////////////// CyFunctionClassCell.proto ////////////////////
static int __Pyx_CyFunction_InitClassCell(PyObject *cyfunctions, PyObject *classobj);/*proto*/

//////////////////// CyFunctionClassCell ////////////////////
//@requires: CythonFunctionShared

static int __Pyx_CyFunction_InitClassCell(PyObject *cyfunctions, PyObject *classobj) {
    Py_ssize_t i, count = __Pyx_PyList_GET_SIZE(cyfunctions);
    #if !CYTHON_ASSUME_SAFE_SIZE
    if (unlikely(count < 0)) return -1;
    #endif

    for (i = 0; i < count; i++) {
        __pyx_CyFunctionObject *m = (__pyx_CyFunctionObject *)
#if CYTHON_ASSUME_SAFE_MACROS && !CYTHON_AVOID_BORROWED_REFS && !CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS
            PyList_GET_ITEM(cyfunctions, i);
#else
            __Pyx_PySequence_ITEM(cyfunctions, i);
        if (unlikely(!m))
            return -1;
#endif
        __Pyx_CyFunction_SetClassObj(m, classobj);
#if !(CYTHON_ASSUME_SAFE_MACROS && !CYTHON_AVOID_BORROWED_REFS && !CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS)
        Py_DECREF((PyObject*)m);
#endif
    }
    return 0;
}

//////////////////// FusedFunction.module_state_decls ////////////////////

PyTypeObject *__pyx_FusedFunctionType;

//////////////////// FusedFunction.module_state_traverse ///////////////////

Py_VISIT(traverse_module_state->__pyx_FusedFunctionType);

//////////////////// FusedFunction.module_state_clear ///////////////////

Py_CLEAR(clear_module_state->__pyx_FusedFunctionType);

//////////////////// FusedFunction.init //////////////////
//@substitute: naming

if (likely(__pyx_FusedFunction_init($module_cname) == 0)); else

//////////////////// FusedFunction.proto ////////////////////

typedef struct {
    __pyx_CyFunctionObject func;
    PyObject *__signatures__;
    PyObject *self;
#if CYTHON_COMPILING_IN_LIMITED_API
    PyMethodDef *ml;
#endif
} __pyx_FusedFunctionObject;

static PyObject *__pyx_FusedFunction_New(PyMethodDef *ml, int flags,
                                         PyObject *qualname, PyObject *closure,
                                         PyObject *module, PyObject *globals,
                                         PyObject *code);

static int __pyx_FusedFunction_clear(__pyx_FusedFunctionObject *self);
static int __pyx_FusedFunction_init(PyObject *module);

#define __Pyx_FusedFunction_USED

//////////////////// FusedFunction ////////////////////
//@requires: CythonFunctionShared
//@substitute: naming

static PyObject *
__pyx_FusedFunction_New(PyMethodDef *ml, int flags,
                        PyObject *qualname, PyObject *closure,
                        PyObject *module, PyObject *globals,
                        PyObject *code)
{
    PyObject *op = __Pyx_CyFunction_Init(
        // __pyx_CyFunctionObject is correct below since that's the cast that we want.
        PyObject_GC_New(__pyx_CyFunctionObject, CGLOBAL(__pyx_FusedFunctionType)),
        ml, flags, qualname, closure, module, globals, code
    );
    if (likely(op)) {
        __pyx_FusedFunctionObject *fusedfunc = (__pyx_FusedFunctionObject *) op;
        fusedfunc->__signatures__ = NULL;
        fusedfunc->self = NULL;
        #if CYTHON_COMPILING_IN_LIMITED_API
        fusedfunc->ml = ml;
        #endif
        PyObject_GC_Track(op);
    }
    return op;
}

static void
__pyx_FusedFunction_dealloc(__pyx_FusedFunctionObject *self)
{
    PyObject_GC_UnTrack(self);
    Py_CLEAR(self->self);
    Py_CLEAR(self->__signatures__);
    __Pyx__CyFunction_dealloc((__pyx_CyFunctionObject *) self);
}

static int
__pyx_FusedFunction_traverse(__pyx_FusedFunctionObject *self,
                             visitproc visit,
                             void *arg)
{
    // Visiting the type is handled in the CyFunction traverse if needed
    Py_VISIT(self->self);
    Py_VISIT(self->__signatures__);
    return __Pyx_CyFunction_traverse((__pyx_CyFunctionObject *) self, visit, arg);
}

static int
__pyx_FusedFunction_clear(__pyx_FusedFunctionObject *self)
{
    Py_CLEAR(self->self);
    Py_CLEAR(self->__signatures__);
    return __Pyx_CyFunction_clear((__pyx_CyFunctionObject *) self);
}


static __pyx_FusedFunctionObject *
__pyx_FusedFunction_descr_get_locked(__pyx_FusedFunctionObject *func, PyObject *obj)
{
    PyObject *module;
    __pyx_FusedFunctionObject *meth;
    #if CYTHON_COMPILING_IN_LIMITED_API
    module = __Pyx_CyFunction_get_module((__pyx_CyFunctionObject *) func, NULL);
    if ((unlikely(!module))) return NULL;
    #else
    module = ((PyCFunctionObject *) func)->m_module;
    #endif

    meth = (__pyx_FusedFunctionObject *) __pyx_FusedFunction_New(
        #if CYTHON_COMPILING_IN_LIMITED_API
                    func->ml,
        #else
                    ((PyCFunctionObject *) func)->m_ml,
        #endif
                    ((__pyx_CyFunctionObject *) func)->flags,
                    ((__pyx_CyFunctionObject *) func)->func_qualname,
                    ((__pyx_CyFunctionObject *) func)->func_closure,
                    module,
                    ((__pyx_CyFunctionObject *) func)->func_globals,
                    ((__pyx_CyFunctionObject *) func)->func_code);
    #if CYTHON_COMPILING_IN_LIMITED_API
    Py_DECREF(module);
    #endif
    if (unlikely(!meth))
        return NULL;

    Py_XINCREF(func->func.defaults);
    meth->func.defaults = func->func.defaults;

    __Pyx_CyFunction_SetClassObj(meth, __Pyx_CyFunction_GetClassObj(func));

    Py_XINCREF(func->__signatures__);
    meth->__signatures__ = func->__signatures__;

    Py_XINCREF(func->func.defaults_tuple);
    meth->func.defaults_tuple = func->func.defaults_tuple;

    Py_XINCREF(obj);
    meth->self = obj;
    return meth;
}

static PyObject *
__pyx_FusedFunction_descr_get(PyObject *self, PyObject *obj, PyObject *type)
{
    __pyx_FusedFunctionObject *func, *meth;

    func = (__pyx_FusedFunctionObject *) self;

    if (func->self || func->func.flags & __Pyx_CYFUNCTION_STATICMETHOD) {
        // Do not allow rebinding and don't do anything for static methods
        Py_INCREF(self);
        return self;
    }

    if (obj == Py_None)
        obj = NULL;

    if (func->func.flags & __Pyx_CYFUNCTION_CLASSMETHOD)
        obj = type;

    if (obj == NULL) {
        // We aren't actually binding to anything, save the effort of rebinding
        Py_INCREF(self);
        return self;
    }

    __Pyx_BEGIN_CRITICAL_SECTION(func);
    meth = __pyx_FusedFunction_descr_get_locked(func, obj);
    __Pyx_END_CRITICAL_SECTION()

    return (PyObject *) meth;
}

static PyObject *
_obj_to_string(PyObject *obj)
{
    if (PyUnicode_CheckExact(obj))
        return __Pyx_NewRef(obj);
    else if (PyType_Check(obj))
        return PyObject_GetAttr(obj, PYIDENT("__name__"));
    else
        return PyObject_Str(obj);
}

static PyObject *
__pyx_FusedFunction_getitem(__pyx_FusedFunctionObject *self, PyObject *idx)
{
    PyObject *signature = NULL;
    PyObject *unbound_result_func;
    PyObject *result_func = NULL;

    if (unlikely(self->__signatures__ == NULL)) {
        PyErr_SetString(PyExc_TypeError, "Function is not fused");
        return NULL;
    }

    if (PyTuple_Check(idx)) {
        Py_ssize_t n = __Pyx_PyTuple_GET_SIZE(idx);
        PyObject *list;
        int i;
        #if !CYTHON_ASSUME_SAFE_SIZE
        if (unlikely(n < 0)) return NULL;
        #endif

        list = PyList_New(n);
        if (unlikely(!list))
            return NULL;

        for (i = 0; i < n; i++) {
            PyObject *string;
#if CYTHON_ASSUME_SAFE_MACROS && !CYTHON_AVOID_BORROWED_REFS
            PyObject *item = PyTuple_GET_ITEM(idx, i);
#else
            PyObject *item = __Pyx_PySequence_ITEM(idx, i);  if (unlikely(!item)) goto __pyx_err;
#endif
            string = _obj_to_string(item);
#if !(CYTHON_ASSUME_SAFE_MACROS && !CYTHON_AVOID_BORROWED_REFS)
            Py_DECREF(item);
#endif
            if (unlikely(!string)) goto __pyx_err;
	    if (__Pyx_PyList_SET_ITEM(list, i, string) < (0)) goto __pyx_err;
        }

        signature = PyUnicode_Join(PYUNICODE("|"), list);
__pyx_err:;
        Py_DECREF(list);
    } else {
        signature = _obj_to_string(idx);
    }

    if (unlikely(!signature))
        return NULL;

    unbound_result_func = PyObject_GetItem(self->__signatures__, signature);

    if (likely(unbound_result_func)) {
        if (self->self) {
            assert(__Pyx_CyFunction_GetClassObj(unbound_result_func) == __Pyx_CyFunction_GetClassObj(self));

            result_func = __pyx_FusedFunction_descr_get(unbound_result_func,
                                                        self->self, self->self);
        } else {
            result_func = unbound_result_func;
            Py_INCREF(result_func);
        }
    }

    Py_DECREF(signature);
    Py_XDECREF(unbound_result_func);

    return result_func;
}

static PyObject *
__pyx_FusedFunction_callfunction(PyObject *func, PyObject *args, PyObject *kw)
{
     __pyx_CyFunctionObject *cyfunc = (__pyx_CyFunctionObject *) func;
    int static_specialized = (cyfunc->flags & __Pyx_CYFUNCTION_STATICMETHOD &&
                              !((__pyx_FusedFunctionObject *) func)->__signatures__);

    if ((cyfunc->flags & __Pyx_CYFUNCTION_CCLASS) && !static_specialized) {
        return __Pyx_CyFunction_CallAsMethod(func, args, kw);
    } else {
        return __Pyx_CyFunction_Call(func, args, kw);
    }
}

// Note: the 'self' from method binding is passed in in the args tuple,
//       whereas PyCFunctionObject's m_self is passed in as the first
//       argument to the C function. For extension methods we need
//       to pass 'self' as 'm_self' and not as the first element of the
//       args tuple.

static PyObject *
__pyx_FusedFunction_call(PyObject *func, PyObject *args, PyObject *kw)
{
    __pyx_FusedFunctionObject *binding_func = (__pyx_FusedFunctionObject *) func;
    Py_ssize_t argc = __Pyx_PyTuple_GET_SIZE(args);
    PyObject *new_args = NULL;
    __pyx_FusedFunctionObject *new_func = NULL;
    PyObject *result = NULL;
    int is_staticmethod = binding_func->func.flags & __Pyx_CYFUNCTION_STATICMETHOD;
    #if !CYTHON_ASSUME_SAFE_SIZE
    if (unlikely(argc < 0)) return NULL;
    #endif

    if (binding_func->self) {
        // Bound method call, put 'self' in the args tuple
        PyObject *self;
        Py_ssize_t i;
        new_args = PyTuple_New(argc + 1);
        if (unlikely(!new_args))
            return NULL;

        self = binding_func->self;

        Py_INCREF(self);
        if (__Pyx_PyTuple_SET_ITEM(new_args, 0, self) < (0)) goto bad;
        self = NULL;

        for (i = 0; i < argc; i++) {
#if CYTHON_ASSUME_SAFE_MACROS && !CYTHON_AVOID_BORROWED_REFS
            PyObject *item = PyTuple_GET_ITEM(args, i);
            Py_INCREF(item);
#else
            PyObject *item = __Pyx_PySequence_ITEM(args, i);  if (unlikely(!item)) goto bad;
#endif
        if (__Pyx_PyTuple_SET_ITEM(new_args, i + 1, item) < (0)) goto bad;
        }

        args = new_args;
    }

    if (binding_func->__signatures__) {
        PyObject *tup;
        if (is_staticmethod && binding_func->func.flags & __Pyx_CYFUNCTION_CCLASS) {
            // FIXME: this seems wrong, but we must currently pass the signatures dict as 'self' argument
            tup = PyTuple_Pack(3, args,
                               kw == NULL ? Py_None : kw,
                               binding_func->func.defaults_tuple);
            if (unlikely(!tup)) goto bad;
            new_func = (__pyx_FusedFunctionObject *) __Pyx_CyFunction_CallMethod(
                func, binding_func->__signatures__, tup, NULL);
        } else {
            tup = PyTuple_Pack(4, binding_func->__signatures__, args,
                               kw == NULL ? Py_None : kw,
                               binding_func->func.defaults_tuple);
            if (unlikely(!tup)) goto bad;
            new_func = (__pyx_FusedFunctionObject *) __pyx_FusedFunction_callfunction(func, tup, NULL);
        }
        Py_DECREF(tup);

        if (unlikely(!new_func))
            goto bad;

        assert(__Pyx_CyFunction_GetClassObj(new_func) == __Pyx_CyFunction_GetClassObj(binding_func));

        func = (PyObject *) new_func;
    }

    result = __pyx_FusedFunction_callfunction(func, args, kw);
bad:
    Py_XDECREF(new_args);
    Py_XDECREF((PyObject *) new_func);
    return result;
}

static PyMemberDef __pyx_FusedFunction_members[] = {
    {"__signatures__",
     T_OBJECT,
     offsetof(__pyx_FusedFunctionObject, __signatures__),
     READONLY,
     0},
    {"__self__", T_OBJECT_EX, offsetof(__pyx_FusedFunctionObject, self), READONLY, 0},
    // For heap-types __module__ appears not to be inherited (so redeclare)
    #if !CYTHON_COMPILING_IN_LIMITED_API
    {"__module__", T_OBJECT, offsetof(PyCFunctionObject, m_module), 0, 0},
    #endif
    {0, 0, 0, 0, 0},
};

static PyGetSetDef __pyx_FusedFunction_getsets[] = {
    // __doc__ is None for the fused function type, but we need it to be
    // a descriptor for the instance's __doc__, so rebuild the descriptor in our subclass
    // (all other descriptors are inherited)
    {"__doc__",  (getter)__Pyx_CyFunction_get_doc, (setter)__Pyx_CyFunction_set_doc, 0, 0},
    // For heap-types __module__ appears not to be inherited (so redeclare)
    #if CYTHON_COMPILING_IN_LIMITED_API
    {"__module__", (getter)__Pyx_CyFunction_get_module, (setter)__Pyx_CyFunction_set_module, 0, 0},
    #endif
    {0, 0, 0, 0, 0}
};

static PyType_Slot __pyx_FusedFunctionType_slots[] = {
    {Py_tp_dealloc, (void *)__pyx_FusedFunction_dealloc},
    {Py_tp_call, (void *)__pyx_FusedFunction_call},
    {Py_tp_traverse, (void *)__pyx_FusedFunction_traverse},
    {Py_tp_clear, (void *)__pyx_FusedFunction_clear},
    {Py_tp_members, (void *)__pyx_FusedFunction_members},
    {Py_tp_getset, (void *)__pyx_FusedFunction_getsets},
    {Py_tp_descr_get, (void *)__pyx_FusedFunction_descr_get},
    {Py_mp_subscript, (void *)__pyx_FusedFunction_getitem},
    {0, 0},
};

static PyType_Spec __pyx_FusedFunctionType_spec = {
    __PYX_TYPE_MODULE_PREFIX "fused_cython_function",
    sizeof(__pyx_FusedFunctionObject),
    0,
    Py_TPFLAGS_IMMUTABLETYPE | Py_TPFLAGS_DISALLOW_INSTANTIATION |
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC , /*tp_flags*/
    __pyx_FusedFunctionType_slots
};

static int __pyx_FusedFunction_init(PyObject *module) {
    $modulestatetype_cname *mstate = __Pyx_PyModule_GetState(module);
    PyObject *bases = PyTuple_Pack(1, mstate->__pyx_CyFunctionType);
    if (unlikely(!bases)) {
        return -1;
    }
    mstate->__pyx_FusedFunctionType = __Pyx_FetchCommonTypeFromSpec(
        mstate->__pyx_CommonTypesMetaclassType, module, &__pyx_FusedFunctionType_spec, bases);
    Py_DECREF(bases);
    if (unlikely(mstate->__pyx_FusedFunctionType == NULL)) {
        return -1;
    }
    return 0;
}

//////////////////// ClassMethod.proto ////////////////////

#if !CYTHON_COMPILING_IN_LIMITED_API
#include "descrobject.h"
#endif
CYTHON_UNUSED static PyObject* __Pyx_Method_ClassMethod(PyObject *method); /*proto*/

//////////////////// ClassMethod ////////////////////
//@requires: ObjectHandling.c::CachedMethodType

static PyObject* __Pyx_Method_ClassMethod(PyObject *method) {
#if CYTHON_COMPILING_IN_PYPY && PYPY_VERSION_NUM <= 0x05080000
    if (PyObject_TypeCheck(method, &PyWrapperDescr_Type)) {
        // cdef classes
        return PyClassMethod_New(method);
    }
#else
#if CYTHON_COMPILING_IN_PYPY
    // special C-API function only in PyPy >= 5.9
    if (PyMethodDescr_Check(method))
#else
    if (__Pyx_TypeCheck(method, &PyMethodDescr_Type))
#endif
    {
#if CYTHON_COMPILING_IN_LIMITED_API
        return PyErr_Format(
            PyExc_SystemError,
            "Cython cannot yet handle classmethod on a MethodDescriptorType (%S) in limited API mode. "
            "This is most likely a classmethod in a cdef class method with binding=False. "
            "Try setting 'binding' to True.",
            method);
#elif CYTHON_COMPILING_IN_GRAAL && defined(GRAALPY_VERSION_NUM) && GRAALPY_VERSION_NUM > 0x19000000
        // cdef classes
        PyTypeObject *d_type = GraalPyDescrObject_GetType(method);
        return PyDescr_NewClassMethod(d_type, GraalPyMethodDescrObject_GetMethod(method));
#elif CYTHON_COMPILING_IN_GRAAL
        // Remove when GraalPy 24 goes EOL
        PyTypeObject *d_type = PyDescrObject_GetType(method);
        return PyDescr_NewClassMethod(d_type, PyMethodDescrObject_GetMethod(method));
#else
        // cdef classes
        PyMethodDescrObject *descr = (PyMethodDescrObject *)method;
        PyTypeObject *d_type = descr->d_common.d_type;
        return PyDescr_NewClassMethod(d_type, descr->d_method);
#endif
    }
#endif
#if !CYTHON_COMPILING_IN_LIMITED_API
    else if (PyMethod_Check(method)) {
        // python classes
        return PyClassMethod_New(PyMethod_GET_FUNCTION(method));
    }
    else {
        return PyClassMethod_New(method);
    }
#else
    {
        PyObject *func=NULL;
        PyObject *builtins, *classmethod, *classmethod_str, *result=NULL;
        if (__Pyx_TypeCheck(method, CGLOBAL(__Pyx_CachedMethodType))) {
            func = PyObject_GetAttrString(method, "__func__");
            if (!func) goto bad;
        } else {
            func = method;
            Py_INCREF(func);
        }
        builtins = PyEval_GetBuiltins(); // borrowed
        if (unlikely(!builtins)) goto bad;
        classmethod_str = PyUnicode_FromString("classmethod");
        if (unlikely(!classmethod_str)) goto bad;
        classmethod = PyObject_GetItem(builtins, classmethod_str);
        Py_DECREF(classmethod_str);
        if (unlikely(!classmethod)) goto bad;
        result = PyObject_CallFunctionObjArgs(classmethod, func, NULL);
        Py_DECREF(classmethod);

        bad:
        Py_XDECREF(func);
        return result;
    }
#endif

}
