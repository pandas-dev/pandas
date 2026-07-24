// This is copied from genobject.c in CPython 3.6.
// Try to keep it in sync by doing this from time to time:
//    sed -e 's|__pyx_||ig'  Cython/Utility/AsyncGen.c | diff -udw - cpython/Objects/genobject.c | less

//////////////////// AsyncGenerator.module_state_decls ////////////////////

// Freelists boost performance 6-10%; they also reduce memory
// fragmentation, as _PyAsyncGenWrappedValue and PyAsyncGenASend
// are short-living objects that are instantiated for every
// __anext__ call.

PyTypeObject *__pyx__PyAsyncGenWrappedValueType;
PyTypeObject *__pyx__PyAsyncGenASendType;
PyTypeObject *__pyx__PyAsyncGenAThrowType;
PyTypeObject *__pyx_AsyncGenType;

#if CYTHON_USE_FREELISTS
struct __pyx__PyAsyncGenWrappedValue *__Pyx_ag_value_freelist[_PyAsyncGen_MAXFREELIST];
int __Pyx_ag_value_freelist_free;

struct __pyx_PyAsyncGenASend *__Pyx_ag_asend_freelist[_PyAsyncGen_MAXFREELIST];
int __Pyx_ag_asend_freelist_free;
#endif

//////////////////// AsyncGenerator.module_state_traverse ////////////////////

Py_VISIT(traverse_module_state->__pyx__PyAsyncGenWrappedValueType);
Py_VISIT(traverse_module_state->__pyx__PyAsyncGenASendType);
Py_VISIT(traverse_module_state->__pyx__PyAsyncGenAThrowType);
Py_VISIT(traverse_module_state->__pyx_AsyncGenType);

//////////////////// AsyncGenerator.module_state_clear ////////////////////

Py_CLEAR(clear_module_state->__pyx__PyAsyncGenWrappedValueType);
Py_CLEAR(clear_module_state->__pyx__PyAsyncGenASendType);
Py_CLEAR(clear_module_state->__pyx__PyAsyncGenAThrowType);
Py_CLEAR(clear_module_state->__pyx_AsyncGenType);

//////////////////// AsyncGenerator.init //////////////////
//@substitute: naming

if (likely(__pyx_AsyncGen_init($module_cname) == 0)); else

//////////////////// AsyncGenerator.proto ////////////////////
//@requires: Coroutine.c::Coroutine

#define __Pyx_AsyncGen_USED
typedef struct {
    __pyx_CoroutineObject coro;
    PyObject *ag_finalizer;
    int ag_hooks_inited;
    int ag_closed;
    int ag_running_async;
} __pyx_PyAsyncGenObject;

#define __Pyx_AsyncGen_CheckExact(obj) __Pyx_IS_TYPE(obj, CGLOBAL(__pyx_AsyncGenType))
#define __pyx_PyAsyncGenASend_CheckExact(o) \
                    __Pyx_IS_TYPE(o, CGLOBAL(__pyx__PyAsyncGenASendType))
#define __pyx_PyAsyncGenAThrow_CheckExact(o) \
                    __Pyx_IS_TYPE(o, CGLOBAL(__pyx__PyAsyncGenAThrowType))

static PyObject *__Pyx_async_gen_anext(PyObject *o);
static CYTHON_INLINE PyObject *__Pyx_async_gen_asend_iternext(PyObject *o);
static PyObject *__Pyx_async_gen_asend_send(PyObject *g, PyObject *arg);
static PyObject *__Pyx_async_gen_asend_close(PyObject *o, PyObject *args);
static PyObject *__Pyx_async_gen_athrow_close(PyObject *o, PyObject *args);

static PyObject *__Pyx__PyAsyncGenValueWrapperNew(PyObject *val);


static __pyx_CoroutineObject *__Pyx_AsyncGen_New(
            __pyx_coroutine_body_t body, PyObject *code, PyObject *closure,
            PyObject *name, PyObject *qualname, PyObject *module_name);

static int __pyx_AsyncGen_init(PyObject *module);
static void __Pyx_PyAsyncGen_Fini(void);

#if CYTHON_USE_FREELISTS && !defined(_PyAsyncGen_MAXFREELIST)
#define _PyAsyncGen_MAXFREELIST 80
#endif

struct __pyx__PyAsyncGenWrappedValue;
struct __pyx_PyAsyncGenASend;

//////////////////// AsyncGenerator.cleanup ////////////////////

__Pyx_PyAsyncGen_Fini();

//////////////////// AsyncGenerator ////////////////////
//@substitute: naming
//@requires: Coroutine.c::Coroutine
//@requires: Coroutine.c::ReturnWithStopIteration
//@requires: ObjectHandling.c::PyObjectCall2Args
//@requires: ExtensionTypes.c::CallTypeTraverse

PyDoc_STRVAR(__Pyx_async_gen_send_doc,
"send(arg) -> send 'arg' into generator,\n\
return next yielded value or raise StopIteration.");

PyDoc_STRVAR(__Pyx_async_gen_close_doc,
"close() -> raise GeneratorExit inside generator.");

PyDoc_STRVAR(__Pyx_async_gen_throw_doc,
"throw(typ[,val[,tb]]) -> raise exception in generator,\n\
return next yielded value or raise StopIteration.");

PyDoc_STRVAR(__Pyx_async_gen_await_doc,
"__await__() -> return a representation that can be passed into the 'await' expression.");

static __pyx_CoroutineObject *__Pyx_AsyncGen_New(
            __pyx_coroutine_body_t body, PyObject *code, PyObject *closure,
            PyObject *name, PyObject *qualname, PyObject *module_name) {
    __pyx_PyAsyncGenObject *gen = PyObject_GC_New(__pyx_PyAsyncGenObject, CGLOBAL(__pyx_AsyncGenType));
    if (unlikely(!gen))
        return NULL;
    gen->ag_finalizer = NULL;
    gen->ag_closed = 0;
    gen->ag_hooks_inited = 0;
    gen->ag_running_async = 0;
    return __Pyx__Coroutine_NewInit((__pyx_CoroutineObject*)gen, body, code, closure, name, qualname, module_name);
}

// COPY STARTS HERE:

static PyObject *__Pyx_async_gen_asend_new(__pyx_PyAsyncGenObject *, PyObject *);
static PyObject *__Pyx_async_gen_athrow_new(__pyx_PyAsyncGenObject *, PyObject *);

static const char *__Pyx_NON_INIT_CORO_MSG = "can't send non-None value to a just-started coroutine";
static const char *__Pyx_ASYNC_GEN_IGNORED_EXIT_MSG = "async generator ignored GeneratorExit";
static const char *__Pyx_ASYNC_GEN_CANNOT_REUSE_SEND_MSG = "cannot reuse already awaited __anext__()/asend()";
static const char *__Pyx_ASYNC_GEN_CANNOT_REUSE_CLOSE_MSG = "cannot reuse already awaited aclose()/athrow()";

typedef enum {
    __PYX_AWAITABLE_STATE_INIT,   /* new awaitable, has not yet been iterated */
    __PYX_AWAITABLE_STATE_ITER,   /* being iterated */
    __PYX_AWAITABLE_STATE_CLOSED, /* closed */
} __pyx_AwaitableState;

typedef struct __pyx_PyAsyncGenASend {
    PyObject_HEAD
    __pyx_PyAsyncGenObject *ags_gen;

    /* Can be NULL, when in the __anext__() mode (equivalent of "asend(None)") */
    PyObject *ags_sendval;

    __pyx_AwaitableState ags_state;
} __pyx_PyAsyncGenASend;


typedef struct {
    PyObject_HEAD
    __pyx_PyAsyncGenObject *agt_gen;

    /* Can be NULL, when in the "aclose()" mode (equivalent of "athrow(GeneratorExit)") */
    PyObject *agt_args;

    __pyx_AwaitableState agt_state;
} __pyx_PyAsyncGenAThrow;


typedef struct __pyx__PyAsyncGenWrappedValue {
    PyObject_HEAD
    PyObject *agw_val;
} __pyx__PyAsyncGenWrappedValue;

#define __pyx__PyAsyncGenWrappedValue_CheckExact(o) \
                    __Pyx_IS_TYPE(o, CGLOBAL(__pyx__PyAsyncGenWrappedValueType))


static int
__Pyx_async_gen_traverse(__pyx_PyAsyncGenObject *gen, visitproc visit, void *arg)
{
    // visiting the type is handled in the base if needed
    Py_VISIT(gen->ag_finalizer);
    return __Pyx_Coroutine_traverse((__pyx_CoroutineObject*)gen, visit, arg);
}


static PyObject *
__Pyx_async_gen_repr(__pyx_CoroutineObject *o)
{
    // avoid NULL pointer dereference for qualname during garbage collection
    return PyUnicode_FromFormat("<async_generator object %S at %p>",
                                o->gi_qualname ? o->gi_qualname : Py_None, o);
}


static int
__Pyx_async_gen_init_hooks_firstiter(__pyx_PyAsyncGenObject *o, PyObject *firstiter)
{
    PyObject *res;
    // at least asyncio stores methods here => optimise the call
#if CYTHON_UNPACK_METHODS
    PyObject *self;
    if (likely(PyMethod_Check(firstiter)) && likely((self = PyMethod_GET_SELF(firstiter)) != NULL)) {
        PyObject *function = PyMethod_GET_FUNCTION(firstiter);
        res = __Pyx_PyObject_Call2Args(function, self, (PyObject*)o);
    } else
#endif
    res = __Pyx_PyObject_CallOneArg(firstiter, (PyObject*)o);

    Py_DECREF(firstiter);

    if (unlikely(res == NULL))
        return 1;

    Py_DECREF(res);
    return 0;
}


static CYTHON_INLINE int
__Pyx_async_gen_init_hooks_done(__pyx_PyAsyncGenObject *o) {
    return o->ag_hooks_inited != 0;
}

static int
__Pyx_async_gen_init_hooks(__pyx_PyAsyncGenObject *o)
{
#if !CYTHON_COMPILING_IN_PYPY && !CYTHON_COMPILING_IN_LIMITED_API
    PyThreadState *tstate;
#endif
    PyObject *finalizer;
    PyObject *firstiter;

    assert (!__Pyx_async_gen_init_hooks_done(o));
    o->ag_hooks_inited = 1;

#if CYTHON_COMPILING_IN_LIMITED_API
    {
        PyObject *hooks_func = PySys_GetObject("get_asyncgen_hooks");
        if (unlikely(!hooks_func)) {
            PyErr_SetString(
                PyExc_AttributeError,
                "Failed to get 'get_asyncgen_hooks' from sys"
            );
            return 1;
        }
        PyObject *async_gen_hooks = PyObject_CallFunctionObjArgs(hooks_func, NULL);
        if (unlikely(!async_gen_hooks)) return 1;
        firstiter = PySequence_GetItem(async_gen_hooks, 0);
        if (unlikely(!firstiter)) {
            Py_DECREF(async_gen_hooks);
            return 1;
        }
        if (firstiter == Py_None) {
            Py_CLEAR(firstiter);
        }

        finalizer = PySequence_GetItem(async_gen_hooks, 1);
        Py_DECREF(async_gen_hooks);

        if (unlikely(!finalizer)) {
            Py_XDECREF(firstiter);
            return 1;
        }
        if (finalizer == Py_None) {
            Py_CLEAR(finalizer);
        }
    }
#endif

#if CYTHON_COMPILING_IN_PYPY
    finalizer = _PyEval_GetAsyncGenFinalizer();
#elif !CYTHON_COMPILING_IN_LIMITED_API
    tstate = __Pyx_PyThreadState_Current;
    finalizer = tstate->async_gen_finalizer;
#endif
    if (finalizer) {
#if !CYTHON_COMPILING_IN_LIMITED_API
        Py_INCREF(finalizer);
#endif
        o->ag_finalizer = finalizer;
    }

#if CYTHON_COMPILING_IN_PYPY
    firstiter = _PyEval_GetAsyncGenFirstiter();
#elif !CYTHON_COMPILING_IN_LIMITED_API
    firstiter = tstate->async_gen_firstiter;
#endif
    if (firstiter) {
#if !CYTHON_COMPILING_IN_LIMITED_API
        Py_INCREF(firstiter);
#endif
        // Transfers the reference.
        if (unlikely(__Pyx_async_gen_init_hooks_firstiter(o, firstiter)))
            return 1;
    }

    return 0;
}


static PyObject *
__Pyx_async_gen_anext(PyObject *g)
{
    __pyx_PyAsyncGenObject *o = (__pyx_PyAsyncGenObject*) g;
    if (!__Pyx_async_gen_init_hooks_done(o) && unlikely(__Pyx_async_gen_init_hooks(o))) {
        return NULL;
    }
    return __Pyx_async_gen_asend_new(o, NULL);
}

static PyObject *
__Pyx_async_gen_anext_method(PyObject *g, PyObject *arg) {
    CYTHON_UNUSED_VAR(arg);
    return __Pyx_async_gen_anext(g);
}


static PyObject *
__Pyx_async_gen_asend(__pyx_PyAsyncGenObject *o, PyObject *arg)
{
    if (!__Pyx_async_gen_init_hooks_done(o) && unlikely(__Pyx_async_gen_init_hooks(o))) {
        return NULL;
    }
    return __Pyx_async_gen_asend_new(o, arg);
}


static PyObject *
__Pyx_async_gen_aclose(__pyx_PyAsyncGenObject *o, PyObject *arg)
{
    CYTHON_UNUSED_VAR(arg);
    if (!__Pyx_async_gen_init_hooks_done(o) && unlikely(__Pyx_async_gen_init_hooks(o))) {
        return NULL;
    }
    return __Pyx_async_gen_athrow_new(o, NULL);
}


static PyObject *
__Pyx_async_gen_athrow(__pyx_PyAsyncGenObject *o, PyObject *args)
{
    if (!__Pyx_async_gen_init_hooks_done(o) && unlikely(__Pyx_async_gen_init_hooks(o))) {
        return NULL;
    }
    return __Pyx_async_gen_athrow_new(o, args);
}


static PyObject *
__Pyx_async_gen_self_method(PyObject *g, PyObject *arg) {
    CYTHON_UNUSED_VAR(arg);
    return __Pyx_NewRef(g);
}


static PyObject *
__Pyx_async_gen_get_running(__pyx_PyAsyncGenObject *o, void *context) {
    CYTHON_UNUSED_VAR(context);
    return __Pyx_NewRef(o->ag_running_async ? Py_True : Py_False);
}


static PyGetSetDef __Pyx_async_gen_getsetlist[] = {
    {"__name__", (getter)__Pyx_Coroutine_get_name, (setter)__Pyx_Coroutine_set_name,
     PyDoc_STR("name of the async generator"), 0},
    {"__qualname__", (getter)__Pyx_Coroutine_get_qualname, (setter)__Pyx_Coroutine_set_qualname,
     PyDoc_STR("qualified name of the async generator"), 0},
    {"ag_running", (getter)__Pyx_async_gen_get_running, NULL, NULL, 0},
    //REMOVED: {(char*) "ag_await", (getter)coro_get_cr_await, NULL,
    //REMOVED:  (char*) PyDoc_STR("object being awaited on, or None")},
    {0, 0, 0, 0, 0} /* Sentinel */
};

static PyMemberDef __Pyx_async_gen_memberlist[] = {
    //REMOVED: {(char*) "ag_frame",   T_OBJECT, offsetof(__pyx_PyAsyncGenObject, ag_frame),   READONLY},
    //REMOVED: {(char*) "ag_code",    T_OBJECT, offsetof(__pyx_PyAsyncGenObject, ag_code),    READONLY},
    //ADDED: "ag_await"
    {"ag_await", T_OBJECT, offsetof(__pyx_CoroutineObject, yieldfrom), READONLY,
     PyDoc_STR("object being awaited on, or None")},
    {"__module__", T_OBJECT, offsetof(__pyx_CoroutineObject, gi_modulename), 0, 0},
#if PY_VERSION_HEX < 0x030C0000 || CYTHON_COMPILING_IN_LIMITED_API
    {"__weaklistoffset__", T_PYSSIZET, offsetof(__pyx_CoroutineObject, gi_weakreflist), READONLY, 0},
#endif
    {0, 0, 0, 0, 0}      /* Sentinel */
};

PyDoc_STRVAR(__Pyx_async_aclose_doc,
"aclose() -> raise GeneratorExit inside generator.");

PyDoc_STRVAR(__Pyx_async_asend_doc,
"asend(v) -> send 'v' in generator.");

PyDoc_STRVAR(__Pyx_async_athrow_doc,
"athrow(typ[,val[,tb]]) -> raise exception in generator.");

PyDoc_STRVAR(__Pyx_async_aiter_doc,
"__aiter__(v) -> return an asynchronous iterator.");

PyDoc_STRVAR(__Pyx_async_anext_doc,
"__anext__(v) -> continue asynchronous iteration and return the next element.");

static PyMethodDef __Pyx_async_gen_methods[] = {
    {"asend", (PyCFunction)__Pyx_async_gen_asend, METH_O, __Pyx_async_asend_doc},
    {"athrow",(PyCFunction)__Pyx_async_gen_athrow, METH_VARARGS, __Pyx_async_athrow_doc},
    {"aclose", (PyCFunction)__Pyx_async_gen_aclose, METH_NOARGS, __Pyx_async_aclose_doc},
    {"__aiter__", (PyCFunction)__Pyx_async_gen_self_method, METH_NOARGS, __Pyx_async_aiter_doc},
    {"__anext__", (PyCFunction)__Pyx_async_gen_anext_method, METH_NOARGS, __Pyx_async_anext_doc},
    {"__reduce_ex__", (PyCFunction) __Pyx_Coroutine_fail_reduce_ex, METH_O, 0},
    {"__reduce__", (PyCFunction) __Pyx_Coroutine_fail_reduce_ex, METH_NOARGS, 0},
    {0, 0, 0, 0}        /* Sentinel */
};


static PyType_Slot __pyx_AsyncGenType_slots[] = {
    {Py_tp_dealloc, (void *)__Pyx_Coroutine_dealloc},
    {Py_am_aiter, (void *)PyObject_SelfIter},
    {Py_am_anext, (void *)__Pyx_async_gen_anext},
    {Py_tp_repr, (void *)__Pyx_async_gen_repr},
    {Py_tp_traverse, (void *)__Pyx_async_gen_traverse},
    {Py_tp_methods, (void *)__Pyx_async_gen_methods},
    {Py_tp_members, (void *)__Pyx_async_gen_memberlist},
    {Py_tp_getset, (void *)__Pyx_async_gen_getsetlist},
#if CYTHON_USE_TP_FINALIZE
    {Py_tp_finalize, (void *)__Pyx_Coroutine_del},
#endif
    {0, 0},
};

static PyType_Spec __pyx_AsyncGenType_spec = {
    __PYX_TYPE_MODULE_PREFIX "async_generator",
    sizeof(__pyx_PyAsyncGenObject),
    0,
#if PY_VERSION_HEX >= 0x030C0000 && !CYTHON_COMPILING_IN_LIMITED_API
    Py_TPFLAGS_MANAGED_WEAKREF |
#endif
    Py_TPFLAGS_IMMUTABLETYPE | Py_TPFLAGS_DISALLOW_INSTANTIATION |
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
    __pyx_AsyncGenType_slots
};


static int
__Pyx_PyAsyncGen_ClearFreeLists(void)
{
    #if CYTHON_USE_FREELISTS
    int ret = CGLOBAL(__Pyx_ag_value_freelist_free) + CGLOBAL(__Pyx_ag_asend_freelist_free);

    while (CGLOBAL(__Pyx_ag_value_freelist_free)) {
        __pyx__PyAsyncGenWrappedValue *o;
        o = CGLOBAL(__Pyx_ag_value_freelist)[--CGLOBAL(__Pyx_ag_value_freelist_free]);
        assert(__pyx__PyAsyncGenWrappedValue_CheckExact((PyObject*)o));
        __Pyx_PyHeapTypeObject_GC_Del(o);
    }

    while (CGLOBAL(__Pyx_ag_asend_freelist_free)) {
        __pyx_PyAsyncGenASend *o;
        o = CGLOBAL(__Pyx_ag_asend_freelist)[--CGLOBAL(__Pyx_ag_asend_freelist_free)];
        assert(__Pyx_IS_TYPE((PyObject*)o, CGLOBAL(__pyx__PyAsyncGenASendType)));
        __Pyx_PyHeapTypeObject_GC_Del(o);
    }

    return ret;
    #else
    return 0;
    #endif
}

static void
__Pyx_PyAsyncGen_Fini(void)
{
    __Pyx_PyAsyncGen_ClearFreeLists();
}


static PyObject *
__Pyx_async_gen_unwrap_value(__pyx_PyAsyncGenObject *gen, PyObject *result, int iternext)
{
    if (result == NULL) {
        PyObject *exc_type = PyErr_Occurred();
        if (!exc_type) {
            PyErr_SetNone(PyExc_StopAsyncIteration);
            gen->ag_closed = 1;
        } else if (__Pyx_PyErr_GivenExceptionMatches2(exc_type, PyExc_StopAsyncIteration, PyExc_GeneratorExit)) {
            gen->ag_closed = 1;
        }

        gen->ag_running_async = 0;
        return NULL;
    }

    if (__pyx__PyAsyncGenWrappedValue_CheckExact(result)) {
        /* async yield */
        __Pyx_ReturnWithStopIteration(((__pyx__PyAsyncGenWrappedValue*)result)->agw_val, 0, iternext);
        Py_DECREF(result);
        gen->ag_running_async = 0;
        return NULL;
    }

    return result;
}


/* ---------- Async Generator ASend Awaitable ------------ */


static void
__Pyx_async_gen_asend_dealloc(__pyx_PyAsyncGenASend *o)
{
    PyObject_GC_UnTrack((PyObject *)o);
    Py_CLEAR(o->ags_gen);
    Py_CLEAR(o->ags_sendval);
    #if CYTHON_USE_FREELISTS
    if (likely(CGLOBAL(__Pyx_ag_asend_freelist_free) < _PyAsyncGen_MAXFREELIST)) {
        assert(__pyx_PyAsyncGenASend_CheckExact((PyObject*)o));
        CGLOBAL(__Pyx_ag_asend_freelist)[CGLOBAL(__Pyx_ag_asend_freelist_free)++] = o;
    } else
    #endif
    {
        __Pyx_PyHeapTypeObject_GC_Del(o);
    }
}

static int
__Pyx_async_gen_asend_traverse(__pyx_PyAsyncGenASend *o, visitproc visit, void *arg)
{
    {
        int e = __Pyx_call_type_traverse((PyObject*)o, 1, visit, arg);
        if (e) return e;
    }
    Py_VISIT(o->ags_gen);
    Py_VISIT(o->ags_sendval);
    return 0;
}

static PyObject *
__Pyx_async_gen_asend_send_impl(PyObject *g, PyObject *arg, int iternext)
{
    __pyx_PyAsyncGenASend *o = (__pyx_PyAsyncGenASend*) g;
    PyObject *retval;

    if (unlikely(o->ags_state == __PYX_AWAITABLE_STATE_CLOSED)) {
        PyErr_SetString(PyExc_RuntimeError, __Pyx_ASYNC_GEN_CANNOT_REUSE_SEND_MSG);
        return NULL;
    }

    if (o->ags_state == __PYX_AWAITABLE_STATE_INIT) {
        if (unlikely(o->ags_gen->ag_running_async)) {
            PyErr_SetString(
                PyExc_RuntimeError,
                "anext(): asynchronous generator is already running");
            return NULL;
        }

        if (arg == NULL || arg == Py_None) {
            arg = o->ags_sendval ? o->ags_sendval : Py_None;
        }
        o->ags_state = __PYX_AWAITABLE_STATE_ITER;
    }

    o->ags_gen->ag_running_async = 1;
    retval = __Pyx_Coroutine_Send((PyObject*)o->ags_gen, arg);
    retval = __Pyx_async_gen_unwrap_value(o->ags_gen, retval, iternext);

    if (!retval) {
        o->ags_state = __PYX_AWAITABLE_STATE_CLOSED;
    }

    return retval;
}

static PyObject *
__Pyx_async_gen_asend_send(PyObject *g, PyObject *arg)
{
    return __Pyx_async_gen_asend_send_impl(g, arg, 0);
}


static CYTHON_INLINE PyObject *
__Pyx_async_gen_asend_iternext(PyObject *o)
{
    return __Pyx_async_gen_asend_send_impl(o, Py_None, 1);
}


static PyObject *
__Pyx_async_gen_asend_throw(__pyx_PyAsyncGenASend *o, PyObject *args)
{
    PyObject *result;

    if (unlikely(o->ags_state == __PYX_AWAITABLE_STATE_CLOSED)) {
        PyErr_SetString(PyExc_RuntimeError, __Pyx_ASYNC_GEN_CANNOT_REUSE_SEND_MSG);
        return NULL;
    }

    result = __Pyx_Coroutine_Throw((PyObject*)o->ags_gen, args);
    result = __Pyx_async_gen_unwrap_value(o->ags_gen, result, 0);

    if (result == NULL) {
        o->ags_state = __PYX_AWAITABLE_STATE_CLOSED;
    }

    return result;
}


static PyObject *
__Pyx_async_gen_asend_close(PyObject *g, PyObject *args)
{
    __pyx_PyAsyncGenASend *o = (__pyx_PyAsyncGenASend*) g;
    CYTHON_UNUSED_VAR(args);
    o->ags_state = __PYX_AWAITABLE_STATE_CLOSED;
    Py_RETURN_NONE;
}


static PyMethodDef __Pyx_async_gen_asend_methods[] = {
    {"send", (PyCFunction)__Pyx_async_gen_asend_send, METH_O, __Pyx_async_gen_send_doc},
    {"throw", (PyCFunction)__Pyx_async_gen_asend_throw, METH_VARARGS, __Pyx_async_gen_throw_doc},
    {"close", (PyCFunction)__Pyx_async_gen_asend_close, METH_NOARGS, __Pyx_async_gen_close_doc},
    {"__await__", (PyCFunction)__Pyx_async_gen_self_method, METH_NOARGS, __Pyx_async_gen_await_doc},
    {0, 0, 0, 0}        /* Sentinel */
};


static PyType_Slot __pyx__PyAsyncGenASendType_slots[] = {
    {Py_tp_dealloc, (void *)__Pyx_async_gen_asend_dealloc},
    {Py_am_await, (void *)PyObject_SelfIter},
    {Py_tp_traverse, (void *)__Pyx_async_gen_asend_traverse},
    {Py_tp_methods, (void *)__Pyx_async_gen_asend_methods},
    {Py_tp_iter, (void *)PyObject_SelfIter},
    {Py_tp_iternext, (void *)__Pyx_async_gen_asend_iternext},
    {0, 0},
};

static PyType_Spec __pyx__PyAsyncGenASendType_spec = {
    __PYX_TYPE_MODULE_PREFIX "async_generator_asend",
    sizeof(__pyx_PyAsyncGenASend),
    0,
    Py_TPFLAGS_IMMUTABLETYPE | Py_TPFLAGS_DISALLOW_INSTANTIATION |
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
    __pyx__PyAsyncGenASendType_slots
};


static PyObject *
__Pyx_async_gen_asend_new(__pyx_PyAsyncGenObject *gen, PyObject *sendval)
{
    __pyx_PyAsyncGenASend *o;
    #if CYTHON_USE_FREELISTS
    if (likely(CGLOBAL(__Pyx_ag_asend_freelist_free))) {
        CGLOBAL(__Pyx_ag_asend_freelist_free)--;
        o = CGLOBAL(__Pyx_ag_asend_freelist)[CGLOBAL(__Pyx_ag_asend_freelist_free)];
        #if CYTHON_COMPILING_IN_LIMITED_API
        Py_DECREF(Py_TYPE((PyObject*)o)); // PyObject_Init resets and increfs the type
        (void) PyObject_Init((PyObject *)o, CGLOBAL(__pyx__PyAsyncGenASendType));
        #else
        _Py_NewReference((PyObject *)o);
        #endif
    } else
    #endif
    {
        o = PyObject_GC_New(__pyx_PyAsyncGenASend, CGLOBAL(__pyx__PyAsyncGenASendType));
        if (unlikely(o == NULL)) {
            return NULL;
        }
    }

    Py_INCREF((PyObject*)gen);
    o->ags_gen = gen;

    Py_XINCREF(sendval);
    o->ags_sendval = sendval;

    o->ags_state = __PYX_AWAITABLE_STATE_INIT;

    PyObject_GC_Track((PyObject*)o);
    return (PyObject*)o;
}


/* ---------- Async Generator Value Wrapper ------------ */


static void
__Pyx_async_gen_wrapped_val_dealloc(__pyx__PyAsyncGenWrappedValue *o)
{
    PyObject_GC_UnTrack((PyObject *)o);
    Py_CLEAR(o->agw_val);
    #if CYTHON_USE_FREELISTS
    if (likely(CGLOBAL(__Pyx_ag_value_freelist_free) < _PyAsyncGen_MAXFREELIST)) {
        assert(__pyx__PyAsyncGenWrappedValue_CheckExact((PyObject*)o));
        CGLOBAL(__Pyx_ag_value_freelist)[CGLOBAL(__Pyx_ag_value_freelist_free)++] = o;
    } else
    #endif
    {
        __Pyx_PyHeapTypeObject_GC_Del(o);
    }
}


static int
__Pyx_async_gen_wrapped_val_traverse(__pyx__PyAsyncGenWrappedValue *o,
                                     visitproc visit, void *arg)
{
    {
        int e = __Pyx_call_type_traverse((PyObject*)o, 1, visit, arg);
        if (e) return e;
    }
    Py_VISIT(o->agw_val);
    return 0;
}


static PyType_Slot __pyx__PyAsyncGenWrappedValueType_slots[] = {
    {Py_tp_dealloc, (void *)__Pyx_async_gen_wrapped_val_dealloc},
    {Py_tp_traverse, (void *)__Pyx_async_gen_wrapped_val_traverse},
    {0, 0},
};

static PyType_Spec __pyx__PyAsyncGenWrappedValueType_spec = {
    __PYX_TYPE_MODULE_PREFIX "async_generator_wrapped_value",
    sizeof(__pyx__PyAsyncGenWrappedValue),
    0,
    Py_TPFLAGS_IMMUTABLETYPE | Py_TPFLAGS_DISALLOW_INSTANTIATION |
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
    __pyx__PyAsyncGenWrappedValueType_slots
};


static PyObject *
__Pyx__PyAsyncGenValueWrapperNew(PyObject *val)
{
    // NOTE: steals a reference to val !
    __pyx__PyAsyncGenWrappedValue *o;
    assert(val);

    #if CYTHON_USE_FREELISTS
    if (likely(CGLOBAL(__Pyx_ag_value_freelist_free))) {
        CGLOBAL(__Pyx_ag_value_freelist_free)--;
        o = CGLOBAL(__Pyx_ag_value_freelist)[CGLOBAL(__Pyx_ag_value_freelist_free)];
        assert(__pyx__PyAsyncGenWrappedValue_CheckExact((PyObject*)o));
        #if CYTHON_COMPILING_IN_LIMITED_API
        Py_DECREF(Py_TYPE((PyObject*)o)); // PyObject_Init resets and increfs the type
        (void) PyObject_Init((PyObject*)o, CGLOBAL(__pyx__PyAsyncGenWrappedValueType));
        #else
        _Py_NewReference((PyObject*)o);
        #endif
    } else
    #endif
    {
        o = PyObject_GC_New(__pyx__PyAsyncGenWrappedValue, CGLOBAL(__pyx__PyAsyncGenWrappedValueType));
        if (unlikely(!o)) {
            Py_DECREF(val);
            return NULL;
        }
    }
    o->agw_val = val;
    // no Py_INCREF(val) - steals reference!
    PyObject_GC_Track((PyObject*)o);
    return (PyObject*)o;
}


/* ---------- Async Generator AThrow awaitable ------------ */


static void
__Pyx_async_gen_athrow_dealloc(__pyx_PyAsyncGenAThrow *o)
{
    PyObject_GC_UnTrack((PyObject *)o);
    Py_CLEAR(o->agt_gen);
    Py_CLEAR(o->agt_args);
    __Pyx_PyHeapTypeObject_GC_Del(o);
}


static int
__Pyx_async_gen_athrow_traverse(__pyx_PyAsyncGenAThrow *o, visitproc visit, void *arg)
{
    {
        int e = __Pyx_call_type_traverse((PyObject*)o, 1, visit, arg);
        if (e) return e;
    }
    Py_VISIT(o->agt_gen);
    Py_VISIT(o->agt_args);
    return 0;
}


static PyObject *
__Pyx_async_gen_athrow_send_impl(__pyx_PyAsyncGenAThrow *o, PyObject *arg, int iternext)
{
    __pyx_CoroutineObject *gen = (__pyx_CoroutineObject*)o->agt_gen;
    PyObject *retval, *exc_type;

    if (unlikely(o->agt_state == __PYX_AWAITABLE_STATE_CLOSED)) {
        PyErr_SetString(PyExc_RuntimeError, __Pyx_ASYNC_GEN_CANNOT_REUSE_CLOSE_MSG);
        return NULL;
    }

    if (unlikely(gen->resume_label == -1)) {
        // already run past the end
        o->agt_state = __PYX_AWAITABLE_STATE_CLOSED;
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }

    if (o->agt_state == __PYX_AWAITABLE_STATE_INIT) {
        if (unlikely(o->agt_gen->ag_running_async)) {
            o->agt_state = __PYX_AWAITABLE_STATE_CLOSED;
            if (o->agt_args == NULL) {
                PyErr_SetString(
                    PyExc_RuntimeError,
                    "aclose(): asynchronous generator is already running");
            } else {
                PyErr_SetString(
                    PyExc_RuntimeError,
                    "athrow(): asynchronous generator is already running");
            }
            return NULL;
        }

        if (unlikely(o->agt_gen->ag_closed)) {
            o->agt_state = __PYX_AWAITABLE_STATE_CLOSED;
            PyErr_SetNone(PyExc_StopAsyncIteration);
            return NULL;
        }

        if (unlikely(arg != Py_None)) {
            PyErr_SetString(PyExc_RuntimeError, __Pyx_NON_INIT_CORO_MSG);
            return NULL;
        }

        o->agt_state = __PYX_AWAITABLE_STATE_ITER;
        o->agt_gen->ag_running_async = 1;

        if (o->agt_args == NULL) {
            /* aclose() mode */
            o->agt_gen->ag_closed = 1;

            retval = __Pyx__Coroutine_Throw((PyObject*)gen,
                /* Do not close generator when PyExc_GeneratorExit is passed */
                PyExc_GeneratorExit, NULL, NULL, NULL, 0);

            if (retval && __pyx__PyAsyncGenWrappedValue_CheckExact(retval)) {
                Py_DECREF(retval);
                goto yield_close;
            }
        } else {
            PyObject *typ;
            PyObject *tb = NULL;
            PyObject *val = NULL;

            if (unlikely(!PyArg_UnpackTuple(o->agt_args, "athrow", 1, 3, &typ, &val, &tb))) {
                return NULL;
            }

            retval = __Pyx__Coroutine_Throw((PyObject*)gen,
                /* Do not close generator when PyExc_GeneratorExit is passed */
                typ, val, tb, o->agt_args, 0);
            retval = __Pyx_async_gen_unwrap_value(o->agt_gen, retval, iternext);
        }
        if (retval == NULL) {
            goto check_error;
        }
        return retval;
    }

    assert (o->agt_state == __PYX_AWAITABLE_STATE_ITER);

    retval = __Pyx_Coroutine_Send((PyObject *)gen, arg);
    if (o->agt_args) {
        return __Pyx_async_gen_unwrap_value(o->agt_gen, retval, iternext);
    } else {
        /* aclose() mode */
        if (retval) {
            if (unlikely(__pyx__PyAsyncGenWrappedValue_CheckExact(retval))) {
                Py_DECREF(retval);
                goto yield_close;
            }
            else {
                return retval;
            }
        }
        else {
            goto check_error;
        }
    }

yield_close:
    o->agt_gen->ag_running_async = 0;
    o->agt_state = __PYX_AWAITABLE_STATE_CLOSED;
    PyErr_SetString(
        PyExc_RuntimeError, __Pyx_ASYNC_GEN_IGNORED_EXIT_MSG);
    return NULL;

check_error:
    o->agt_gen->ag_running_async = 0;
    o->agt_state = __PYX_AWAITABLE_STATE_CLOSED;
    exc_type = PyErr_Occurred();
    if (__Pyx_PyErr_GivenExceptionMatches2(exc_type, PyExc_StopAsyncIteration, PyExc_GeneratorExit)) {
        if (o->agt_args == NULL) {
            // when aclose() is called we don't want to propagate
            // StopAsyncIteration or GeneratorExit; just raise
            // StopIteration, signalling that this 'aclose()' await
            // is done.
            PyErr_Clear();
            PyErr_SetNone(PyExc_StopIteration);
        }
    }
    return NULL;
}

static PyObject *
__Pyx_async_gen_athrow_send(__pyx_PyAsyncGenAThrow *o, PyObject *arg)
{
    return __Pyx_async_gen_athrow_send_impl(o, arg, 0);
}


static PyObject *
__Pyx_async_gen_athrow_throw(__pyx_PyAsyncGenAThrow *o, PyObject *args)
{
    PyObject *retval;

    if (unlikely(o->agt_state == __PYX_AWAITABLE_STATE_CLOSED)) {
        PyErr_SetString(PyExc_RuntimeError, __Pyx_ASYNC_GEN_CANNOT_REUSE_CLOSE_MSG);
        return NULL;
    }

    retval = __Pyx_Coroutine_Throw((PyObject*)o->agt_gen, args);
    if (o->agt_args) {
        return __Pyx_async_gen_unwrap_value(o->agt_gen, retval, 0);
    } else {
        // aclose() mode
        PyObject *exc_type;
        if (unlikely(retval && __pyx__PyAsyncGenWrappedValue_CheckExact(retval))) {
            o->agt_gen->ag_running_async = 0;
            o->agt_state = __PYX_AWAITABLE_STATE_CLOSED;
            Py_DECREF(retval);
            PyErr_SetString(PyExc_RuntimeError, __Pyx_ASYNC_GEN_IGNORED_EXIT_MSG);
            return NULL;
        }
        exc_type = PyErr_Occurred();
        if (__Pyx_PyErr_GivenExceptionMatches2(exc_type, PyExc_StopAsyncIteration, PyExc_GeneratorExit)) {
            // when aclose() is called we don't want to propagate
            // StopAsyncIteration or GeneratorExit; just raise
            // StopIteration, signalling that this 'aclose()' await
            // is done.
            PyErr_Clear();
            PyErr_SetNone(PyExc_StopIteration);
        }
        return retval;
    }
}


static PyObject *
__Pyx_async_gen_athrow_iternext(__pyx_PyAsyncGenAThrow *o)
{
    return __Pyx_async_gen_athrow_send_impl(o, Py_None, 1);
}


static PyObject *
__Pyx_async_gen_athrow_close(PyObject *g, PyObject *args)
{
    __pyx_PyAsyncGenAThrow *o = (__pyx_PyAsyncGenAThrow*) g;
    CYTHON_UNUSED_VAR(args);
    o->agt_state = __PYX_AWAITABLE_STATE_CLOSED;
    Py_RETURN_NONE;
}


static PyMethodDef __Pyx_async_gen_athrow_methods[] = {
    {"send", (PyCFunction)__Pyx_async_gen_athrow_send, METH_O, __Pyx_async_gen_send_doc},
    {"throw", (PyCFunction)__Pyx_async_gen_athrow_throw, METH_VARARGS, __Pyx_async_gen_throw_doc},
    {"close", (PyCFunction)__Pyx_async_gen_athrow_close, METH_NOARGS, __Pyx_async_gen_close_doc},
    {"__await__", (PyCFunction)__Pyx_async_gen_self_method, METH_NOARGS, __Pyx_async_gen_await_doc},
    {0, 0, 0, 0}        /* Sentinel */
};


static PyType_Slot __pyx__PyAsyncGenAThrowType_slots[] = {
    {Py_tp_dealloc, (void *)__Pyx_async_gen_athrow_dealloc},
    {Py_am_await, (void *)PyObject_SelfIter},
    {Py_tp_traverse, (void *)__Pyx_async_gen_athrow_traverse},
    {Py_tp_iter, (void *)PyObject_SelfIter},
    {Py_tp_iternext, (void *)__Pyx_async_gen_athrow_iternext},
    {Py_tp_methods, (void *)__Pyx_async_gen_athrow_methods},
    {Py_tp_getattro, (void *)PyObject_GenericGetAttr},
    {0, 0},
};

static PyType_Spec __pyx__PyAsyncGenAThrowType_spec = {
    __PYX_TYPE_MODULE_PREFIX "async_generator_athrow",
    sizeof(__pyx_PyAsyncGenAThrow),
    0,
    Py_TPFLAGS_IMMUTABLETYPE | Py_TPFLAGS_DISALLOW_INSTANTIATION |
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
    __pyx__PyAsyncGenAThrowType_slots
};


static PyObject *
__Pyx_async_gen_athrow_new(__pyx_PyAsyncGenObject *gen, PyObject *args)
{
    __pyx_PyAsyncGenAThrow *o;
    o = PyObject_GC_New(__pyx_PyAsyncGenAThrow, CGLOBAL(__pyx__PyAsyncGenAThrowType));
    if (unlikely(o == NULL)) {
        return NULL;
    }
    o->agt_gen = gen;
    o->agt_args = args;
    o->agt_state = __PYX_AWAITABLE_STATE_INIT;
    Py_INCREF((PyObject*)gen);
    Py_XINCREF(args);
    PyObject_GC_Track((PyObject*)o);
    return (PyObject*)o;
}


/* ---------- global type sharing ------------ */

static int __pyx_AsyncGen_init(PyObject *module) {
    $modulestatetype_cname *mstate = __Pyx_PyModule_GetState(module);
    mstate->__pyx_AsyncGenType = __Pyx_FetchCommonTypeFromSpec(
        mstate->__pyx_CommonTypesMetaclassType, module, &__pyx_AsyncGenType_spec, NULL);
    if (unlikely(!mstate->__pyx_AsyncGenType))
        return -1;

    mstate->__pyx__PyAsyncGenAThrowType = __Pyx_FetchCommonTypeFromSpec(NULL, module, &__pyx__PyAsyncGenAThrowType_spec, NULL);
    if (unlikely(!mstate->__pyx__PyAsyncGenAThrowType))
        return -1;

    mstate->__pyx__PyAsyncGenWrappedValueType = __Pyx_FetchCommonTypeFromSpec(NULL, module, &__pyx__PyAsyncGenWrappedValueType_spec, NULL);
    if (unlikely(!mstate->__pyx__PyAsyncGenWrappedValueType))
        return -1;

    mstate->__pyx__PyAsyncGenASendType = __Pyx_FetchCommonTypeFromSpec(NULL, module, &__pyx__PyAsyncGenASendType_spec, NULL);
    if (unlikely(!mstate->__pyx__PyAsyncGenASendType))
        return -1;

    return 0;
}
