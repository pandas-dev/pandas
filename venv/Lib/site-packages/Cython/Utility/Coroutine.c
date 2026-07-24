// ### CPython revisions to investigate for potential changes:

// 2f180ce2cb6 (bpo-44530: Add co_qualname field to PyCodeObject (GH-26941))
// a4760cc32d9 (bpo-42747: Remove Py_TPFLAGS_HAVE_AM_SEND and make Py_TPFLAGS_HAVE_VERSION_TAG no-op (GH-27260))
// ae0a2b75625 (bpo-44590: Lazily allocate frame objects (GH-27077))
// 69806b9516d (bpo-46009: Do not exhaust generator when send() method raises (GH-29986))
// b04dfbbe4bd (bpo-46409: Make generators in bytecode (GH-30633) -- ".gi_suspended")
// 304197b3820 (bpo-46944: use FASTCALL calling convention in generator.throw (GH-31723))
// 8be7c2bc5ad (bpo-14911: Corrected generator.throw() documentation (GH-32207))
// 5bfb3c372bd (GH-90997: Wrap yield from/await in a virtual try/except StopIteration (GH-96010))
// 83a3de4e063 (gh-96348: Deprecate the 3-arg signature of coroutine.throw and generator.throw (GH-96428))
// f4adb975061 (GH-96793: Implement PEP 479 in bytecode. (GH-99006))
// f02fa64bf2d (GH-100762: Don't call `gen.throw()` in `gen.close()`, unless necessary. (GH-101013))
// 11a2c6ce516 (gh-102192: Replace PyErr_Fetch/Restore etc by more efficient alternatives (in Objects/) (#102218))
// ced13c96a4e (gh-79940: add introspection API for asynchronous generators to `inspect` module (#11590) -- ".ag_suspended")
// d56c933992c (gh-104770: Let generator.close() return value (#104771))
// 7fc542c88dc (GH-89091: raise `RuntimeWarning` for unawaited async generator methods (#104611))
// 0d30a5a4096 (GH-100964: Break cycles involving exception state when returning from generator (GH-107563))
// 7d369d471cf (GH-117536: GH-117894: fix athrow().throw(...) unawaited warning (GH-117851))
// fc7e1aa3c00 (GH-117881: fix athrow().throw()/asend().throw() concurrent access (GH-117882))
// e5c699280de (GH-117714: implement athrow().close() and asend().close() using throw (GH-117906))
// 2c7209a3bdf (gh-114091: Reword error message for unawaitable types (#114090))


//////////////////// CoroutineSetYieldFrom ////////////////////

static void
__Pyx_Coroutine_Set_Owned_Yield_From(__pyx_CoroutineObject *gen, PyObject *yf) {
    // NOTE: steals a reference to yf by transferring it to 'gen->yieldfrom' !
    assert (!gen->yieldfrom);
#if CYTHON_USE_AM_SEND
    assert (!gen->yieldfrom_am_send);
    #if PY_VERSION_HEX < 0x030A00F0
    if (__Pyx_PyType_HasFeature(Py_TYPE(yf), __Pyx_TPFLAGS_HAVE_AM_SEND))
    #endif
    {
        __Pyx_pyiter_sendfunc am_send;
        #if __PYX_LIMITED_VERSION_HEX >= 0x030A0000
        am_send = __Pyx_PyObject_TryGetSubSlot(yf, tp_as_async, am_send, __Pyx_pyiter_sendfunc);
        #else
        __Pyx_PyAsyncMethodsStruct* tp_as_async = (__Pyx_PyAsyncMethodsStruct*) Py_TYPE(yf)->tp_as_async;
        am_send = tp_as_async ? tp_as_async->am_send : NULL;
        #endif
        if (likely(am_send)) {
            // Keep a direct reference to am->am_send, if provided.
            gen->yieldfrom_am_send = am_send;
        }
    }
#endif
    gen->yieldfrom = yf;
}


//////////////////// GeneratorYieldFrom.proto ////////////////////

static CYTHON_INLINE __Pyx_PySendResult __Pyx_Generator_Yield_From(__pyx_CoroutineObject *gen, PyObject *source, PyObject **retval);

//////////////////// GeneratorYieldFrom ////////////////////
//@requires: Generator
//@requires: CoroutineSetYieldFrom
//@requires: ObjectHandling.c::IterNextPlain

#if CYTHON_USE_TYPE_SLOTS
static void __Pyx_PyIter_CheckErrorAndDecref(PyObject *source) {
    __Pyx_TypeName source_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(source));
    PyErr_Format(PyExc_TypeError,
        "iter() returned non-iterator of type '" __Pyx_FMT_TYPENAME "'", source_type_name);
    __Pyx_DECREF_TypeName(source_type_name);
    Py_DECREF(source);
}
#endif

static CYTHON_INLINE __Pyx_PySendResult __Pyx_Generator_Yield_From(__pyx_CoroutineObject *gen, PyObject *source, PyObject **retval) {
    PyObject *source_gen;
    __Pyx_PySendResult result;
#ifdef __Pyx_Coroutine_USED
    if (__Pyx_Coroutine_Check(source)) {
        // TODO: this should only happen for types.coroutine()ed generators, but we can't determine that here
        Py_INCREF(source);
        source_gen = source;
        result = __Pyx_Coroutine_AmSend(source, Py_None, retval);
    } else
#endif
    {
#if CYTHON_USE_TYPE_SLOTS
        if (likely(Py_TYPE(source)->tp_iter)) {
            source_gen = Py_TYPE(source)->tp_iter(source);
            if (unlikely(!source_gen)) {
                *retval = NULL;
                return PYGEN_ERROR;
            }
            if (unlikely(!PyIter_Check(source_gen))) {
                __Pyx_PyIter_CheckErrorAndDecref(source_gen);
                *retval = NULL;
                return PYGEN_ERROR;
            }
        } else
        // CPython also allows non-iterable sequences to be iterated over
#endif
        {
            source_gen = PyObject_GetIter(source);
            if (unlikely(!source_gen)) {
                *retval = NULL;
                return PYGEN_ERROR;
            }
        }
        // source_gen is now the iterator, make the first next() call
        *retval = __Pyx_PyIter_Next_Plain(source_gen);
        result = __Pyx_Coroutine_status_from_result(retval);
    }
    if (likely(result == PYGEN_NEXT)) {
        __Pyx_Coroutine_Set_Owned_Yield_From(gen, source_gen);
        return PYGEN_NEXT;
    }
    Py_DECREF(source_gen);
    return result;
}


//////////////////// CoroutineYieldFrom.proto ////////////////////

static CYTHON_INLINE __Pyx_PySendResult __Pyx_Coroutine_Yield_From(__pyx_CoroutineObject *gen, PyObject *source, PyObject **retval); /*proto*/

//////////////////// CoroutineYieldFrom ////////////////////
//@requires: Coroutine
//@requires: GetAwaitIter
//@requires: CoroutineSetYieldFrom
//@requires: ObjectHandling.c::IterNextPlain

static __Pyx_PySendResult __Pyx_Coroutine_Yield_From_Coroutine(__pyx_CoroutineObject *gen, PyObject *source, PyObject **retval) {
    __Pyx_PySendResult result;
    if (unlikely(((__pyx_CoroutineObject*)source)->yieldfrom)) {
        PyErr_SetString(
            PyExc_RuntimeError,
            "coroutine is being awaited already");
        return PYGEN_ERROR;
    }
    result = __Pyx_Coroutine_AmSend(source, Py_None, retval);
    if (result == PYGEN_NEXT) {
        Py_INCREF(source);
        __Pyx_Coroutine_Set_Owned_Yield_From(gen, source);
    }
    return result;
}

static __Pyx_PySendResult __Pyx_Coroutine_Yield_From_Generic(__pyx_CoroutineObject *gen, PyObject *source, PyObject **retval) {
    __Pyx_PySendResult result;
    PyObject *source_gen = NULL;

    source_gen = __Pyx_Coroutine_GetAwaitableIter(source);
    if (unlikely(!source_gen)) return PYGEN_ERROR;

    // source_gen is now the iterator, make the first next() call
    if (__Pyx_Coroutine_Check(source_gen)) {
        result = __Pyx_Coroutine_Yield_From_Coroutine(gen, source_gen, retval);
        Py_DECREF(source_gen);
        return result;
    }

    *retval = __Pyx_PyIter_Next_Plain(source_gen);
    if (*retval) {
        __Pyx_Coroutine_Set_Owned_Yield_From(gen, source_gen);
        return PYGEN_NEXT;
    }

    result = __Pyx_Coroutine_status_from_result(retval);
    Py_XDECREF(source_gen);
    return result;
}

static CYTHON_INLINE __Pyx_PySendResult __Pyx_Coroutine_Yield_From(__pyx_CoroutineObject *gen, PyObject *source, PyObject **retval) {
    if (__Pyx_Coroutine_Check(source)) {
        return __Pyx_Coroutine_Yield_From_Coroutine(gen, source, retval);
    }

#ifdef __Pyx_AsyncGen_USED
    // Inlined "__pyx_PyAsyncGenASend" handling to avoid the series of generic calls.
    if (__pyx_PyAsyncGenASend_CheckExact(source)) {
        *retval = __Pyx_async_gen_asend_iternext(source);
        if (*retval) {
            Py_INCREF(source);
            __Pyx_Coroutine_Set_Owned_Yield_From(gen, source);
            return PYGEN_NEXT;
        }
        return __Pyx_Coroutine_status_from_result(retval);
    }
#endif

    return __Pyx_Coroutine_Yield_From_Generic(gen, source, retval);
}


//////////////////// GetAwaitIter.proto ////////////////////

static CYTHON_INLINE PyObject *__Pyx_Coroutine_GetAwaitableIter(PyObject *o); /*proto*/
static PyObject *__Pyx__Coroutine_GetAwaitableIter(PyObject *o); /*proto*/

//////////////////// GetAwaitIter ////////////////////
//@requires: ObjectHandling.c::PyObjectGetMethod
//@requires: ObjectHandling.c::PyObjectCallNoArg
//@requires: ObjectHandling.c::PyObjectCallOneArg
//@requires: Coro_CheckExact

static CYTHON_INLINE PyObject *__Pyx_Coroutine_GetAwaitableIter(PyObject *o) {
#ifdef __Pyx_Coroutine_USED
    if (__Pyx_Coroutine_Check(o)) {
        return __Pyx_NewRef(o);
    }
#endif
    return __Pyx__Coroutine_GetAwaitableIter(o);
}


static void __Pyx_Coroutine_AwaitableIterError(PyObject *source) {
#if (PY_VERSION_HEX < 0x030d0000 || defined(_PyErr_FormatFromCause)) && !CYTHON_COMPILING_IN_LIMITED_API
    __Pyx_TypeName source_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(source));
    _PyErr_FormatFromCause(PyExc_TypeError,
        "'async for' received an invalid object from __anext__: " __Pyx_FMT_TYPENAME, source_type_name);
    __Pyx_DECREF_TypeName(source_type_name);
#else
    PyObject *exc, *val, *val2, *tb;
    __Pyx_TypeName source_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(source));
    assert(PyErr_Occurred());
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_NormalizeException(&exc, &val, &tb);
    if (tb != NULL) {
        PyException_SetTraceback(val, tb);
        Py_DECREF(tb);
    }
    Py_DECREF(exc);
    assert(!PyErr_Occurred());
    PyErr_Format(PyExc_TypeError,
        "'async for' received an invalid object from __anext__: " __Pyx_FMT_TYPENAME, source_type_name);
    __Pyx_DECREF_TypeName(source_type_name);

    PyErr_Fetch(&exc, &val2, &tb);
    PyErr_NormalizeException(&exc, &val2, &tb);
    Py_INCREF(val);
    PyException_SetCause(val2, val);
    PyException_SetContext(val2, val);
    PyErr_Restore(exc, val2, tb);
#endif
}

// adapted from genobject.c in Py3.5
static PyObject *__Pyx__Coroutine_GetAwaitableIter(PyObject *obj) {
    PyObject *res;
    unaryfunc am_await;
    am_await = __Pyx_PyObject_TryGetSubSlot(obj, tp_as_async, am_await, unaryfunc);
    if (likely(am_await)) {
        res = (*am_await)(obj);
    } else
#if CYTHON_COMPILING_IN_CPYTHON && defined(CO_ITERABLE_COROUTINE)
#if PY_VERSION_HEX >= 0x030C00A6
    if (PyGen_CheckExact(obj) && (PyGen_GetCode((PyGenObject*)obj)->co_flags & CO_ITERABLE_COROUTINE)) {
#else
    if (PyGen_CheckExact(obj) && ((PyGenObject*)obj)->gi_code && ((PyCodeObject *)((PyGenObject*)obj)->gi_code)->co_flags & CO_ITERABLE_COROUTINE) {
#endif
        // Python generator marked with "@types.coroutine" decorator
        return __Pyx_NewRef(obj);
    } else
#endif
    {
#if CYTHON_VECTORCALL && (__PYX_LIMITED_VERSION_HEX >= 0x030C0000 || (!CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX >= 0x03090000))
        PyObject *args[] = { obj };
        res = PyObject_VectorcallMethod(PYIDENT("__await__"), args, 1 |  PY_VECTORCALL_ARGUMENTS_OFFSET, NULL);
        if (unlikely(!res)) {
            // We need to distinguish between "doesn't have __await__" vs "__await__ failed"
            if (PyErr_ExceptionMatches(PyExc_AttributeError)) {
                PyObject *type, *value, *traceback;
                PyErr_Fetch(&type, &value, &traceback);
                int has_attr = PyObject_HasAttr(obj, PYIDENT("__await__"));
                PyErr_Restore(type, value, traceback);
                if (!has_attr) goto slot_error;
            }
        }
#else
        PyObject *method = NULL;
        int is_method = __Pyx_PyObject_GetMethod(obj, PYIDENT("__await__"), &method);
        if (likely(is_method)) {
            res = __Pyx_PyObject_CallOneArg(method, obj);
        } else if (likely(method)) {
            res = __Pyx_PyObject_CallNoArg(method);
        } else
            goto slot_error;
        Py_DECREF(method);
#endif
    }
    if (unlikely(!res)) {
        // surprisingly, CPython replaces the exception here...
        __Pyx_Coroutine_AwaitableIterError(obj);
        goto bad;
    }
    if (unlikely(!PyIter_Check(res))) {
        __Pyx_TypeName res_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(res));
        PyErr_Format(PyExc_TypeError,
            "__await__() returned non-iterator of type '" __Pyx_FMT_TYPENAME "'", res_type_name);
        __Pyx_DECREF_TypeName(res_type_name);
        Py_CLEAR(res);
    } else {
        int is_coroutine = 0;
        #ifdef __Pyx_Coroutine_USED
        is_coroutine |= __Pyx_Coroutine_Check(res);
        #endif
        is_coroutine |= __Pyx_PyCoro_CheckExact(res);
        if (unlikely(is_coroutine)) {
            /* __await__ must return an *iterator*, not
               a coroutine or another awaitable (see PEP 492) */
            PyErr_SetString(PyExc_TypeError,
                            "__await__() returned a coroutine");
            Py_CLEAR(res);
        }
    }
    return res;
slot_error:
    {
        __Pyx_TypeName obj_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(obj));
        PyErr_Format(PyExc_TypeError,
            "object " __Pyx_FMT_TYPENAME " can't be used in 'await' expression", obj_type_name);
        __Pyx_DECREF_TypeName(obj_type_name);
    }
bad:
    return NULL;
}


//////////////////// AsyncIter.proto ////////////////////

static CYTHON_INLINE PyObject *__Pyx_Coroutine_GetAsyncIter(PyObject *o); /*proto*/
static CYTHON_INLINE PyObject *__Pyx_Coroutine_AsyncIterNext(PyObject *o); /*proto*/

//////////////////// AsyncIter ////////////////////
//@requires: GetAwaitIter

static PyObject *__Pyx_Coroutine_GetAsyncIter_Fail(PyObject *obj) {
    __Pyx_TypeName obj_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(obj));
    PyErr_Format(PyExc_TypeError,
                 "'async for' requires an object with __aiter__ method, got " __Pyx_FMT_TYPENAME, obj_type_name);
    __Pyx_DECREF_TypeName(obj_type_name);
    return NULL;
}


static CYTHON_INLINE PyObject *__Pyx_Coroutine_GetAsyncIter(PyObject *obj) {
#ifdef __Pyx_AsyncGen_USED
    if (__Pyx_AsyncGen_CheckExact(obj)) {
        return __Pyx_NewRef(obj);
    }
#endif
    {
        unaryfunc am_aiter = __Pyx_PyObject_TryGetSubSlot(obj, tp_as_async, am_aiter, unaryfunc);
        if (likely(am_aiter)) {
            return (*am_aiter)(obj);
        }
    }
    return __Pyx_Coroutine_GetAsyncIter_Fail(obj);
}


static PyObject *__Pyx_Coroutine_AsyncIterNext_Fail(PyObject *obj) {
    __Pyx_TypeName obj_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(obj));
    PyErr_Format(PyExc_TypeError,
        "'async for' requires an object with __anext__ method, got " __Pyx_FMT_TYPENAME, obj_type_name);
    __Pyx_DECREF_TypeName(obj_type_name);
    return NULL;
}


static CYTHON_INLINE PyObject *__Pyx_Coroutine_AsyncIterNext(PyObject *obj) {
#ifdef __Pyx_AsyncGen_USED
    if (__Pyx_AsyncGen_CheckExact(obj)) {
        return __Pyx_async_gen_anext(obj);
    }
#endif
    {
        unaryfunc am_anext = __Pyx_PyObject_TryGetSubSlot(obj, tp_as_async, am_anext, unaryfunc);
        if (likely(am_anext)) {
            return (*am_anext)(obj);
        }
    }
    return __Pyx_Coroutine_AsyncIterNext_Fail(obj);
}


//////////////////// pep479.proto ////////////////////

static void __Pyx_Generator_Replace_StopIteration(int in_async_gen); /*proto*/

//////////////////// pep479 ////////////////////
//@requires: Exceptions.c::GetException

static void __Pyx_Generator_Replace_StopIteration(int in_async_gen) {
    PyObject *exc, *val, *tb, *cur_exc, *new_exc;
    __Pyx_PyThreadState_declare
    int is_async_stopiteration = 0;
    CYTHON_MAYBE_UNUSED_VAR(in_async_gen);

    __Pyx_PyThreadState_assign
    cur_exc = __Pyx_PyErr_CurrentExceptionType();
    if (likely(!__Pyx_PyErr_GivenExceptionMatches(cur_exc, PyExc_StopIteration))) {
        if (in_async_gen && unlikely(__Pyx_PyErr_GivenExceptionMatches(cur_exc, PyExc_StopAsyncIteration))) {
            is_async_stopiteration = 1;
        } else {
            return;
        }
    }

    // Explicitly chain the Stop(Async)Iteration by moving it to exc_info before creating the RuntimeError.
    __Pyx_GetException(&exc, &val, &tb);
    Py_XDECREF(exc);
    Py_XDECREF(tb);
    new_exc = PyObject_CallFunction(PyExc_RuntimeError, "s",
        is_async_stopiteration ? "async generator raised StopAsyncIteration" :
        in_async_gen ? "async generator raised StopIteration" :
        "generator raised StopIteration");
    if (!new_exc) {
        Py_XDECREF(val);
        return;
    }
    PyException_SetCause(new_exc, val); // steals ref to val
    PyErr_SetObject(PyExc_RuntimeError, new_exc);
    Py_DECREF(new_exc);
}


//////////////////// CoroutineBase.proto ////////////////////
//@substitute: naming

struct __pyx_CoroutineObject;
typedef PyObject *(*__pyx_coroutine_body_t)(struct __pyx_CoroutineObject *, PyThreadState *, PyObject *);

#if CYTHON_USE_EXC_INFO_STACK
// See  https://bugs.python.org/issue25612
#define __Pyx_ExcInfoStruct  _PyErr_StackItem
#else
// Minimal replacement struct for Py<3.7, without the Py3.7 exception state stack.
typedef struct {
    PyObject *exc_type;
    PyObject *exc_value;
    PyObject *exc_traceback;
} __Pyx_ExcInfoStruct;
#endif

typedef struct __pyx_CoroutineObject {
    PyObject_HEAD
    __pyx_coroutine_body_t body;
    PyObject *closure;
    __Pyx_ExcInfoStruct gi_exc_state;
#if PY_VERSION_HEX < 0x030C0000 || CYTHON_COMPILING_IN_LIMITED_API
    PyObject *gi_weakreflist;
#endif
    PyObject *classobj;
    PyObject *yieldfrom;
    // "yieldfrom_am_send" is always included to avoid changing the struct layout.
    __Pyx_pyiter_sendfunc yieldfrom_am_send;
    PyObject *gi_name;
    PyObject *gi_qualname;
    PyObject *gi_modulename;
    PyObject *gi_code;
    PyObject *gi_frame;
#if CYTHON_USE_SYS_MONITORING && (CYTHON_PROFILE || CYTHON_TRACE)
    PyMonitoringState $monitoring_states_cname[__Pyx_MonitoringEventTypes_CyGen_count];
    uint64_t $monitoring_version_cname;
#endif
    int resume_label;
    // is_running is the main thread-safety mechanism for coroutines, so treat it carefully
    char is_running;
} __pyx_CoroutineObject;

static __pyx_CoroutineObject *__Pyx__Coroutine_New(
    PyTypeObject *type, __pyx_coroutine_body_t body, PyObject *code, PyObject *closure,
    PyObject *name, PyObject *qualname, PyObject *module_name); /*proto*/

static __pyx_CoroutineObject *__Pyx__Coroutine_NewInit(
            __pyx_CoroutineObject *gen, __pyx_coroutine_body_t body, PyObject *code, PyObject *closure,
            PyObject *name, PyObject *qualname, PyObject *module_name); /*proto*/

static CYTHON_INLINE void __Pyx_Coroutine_ExceptionClear(__Pyx_ExcInfoStruct *self);
static int __Pyx_Coroutine_clear(PyObject *self); /*proto*/
static __Pyx_PySendResult __Pyx_Coroutine_AmSend(PyObject *self, PyObject *value, PyObject **retval); /*proto*/
static PyObject *__Pyx_Coroutine_Send(PyObject *self, PyObject *value); /*proto*/
static __Pyx_PySendResult __Pyx_Coroutine_Close(PyObject *self, PyObject **retval); /*proto*/
static PyObject *__Pyx_Coroutine_Throw(PyObject *gen, PyObject *args); /*proto*/

// macros for exception state swapping instead of inline functions to make use of the local thread state context
#if CYTHON_USE_EXC_INFO_STACK
#define __Pyx_Coroutine_SwapException(self)
#define __Pyx_Coroutine_ResetAndClearException(self)  __Pyx_Coroutine_ExceptionClear(&(self)->gi_exc_state)
#else
#define __Pyx_Coroutine_SwapException(self) { \
    __Pyx_ExceptionSwap(&(self)->gi_exc_state.exc_type, &(self)->gi_exc_state.exc_value, &(self)->gi_exc_state.exc_traceback); \
    __Pyx_Coroutine_ResetFrameBackpointer(&(self)->gi_exc_state); \
    }
#define __Pyx_Coroutine_ResetAndClearException(self) { \
    __Pyx_ExceptionReset((self)->gi_exc_state.exc_type, (self)->gi_exc_state.exc_value, (self)->gi_exc_state.exc_traceback); \
    (self)->gi_exc_state.exc_type = (self)->gi_exc_state.exc_value = (self)->gi_exc_state.exc_traceback = NULL; \
    }
#endif

#if CYTHON_FAST_THREAD_STATE
#define __Pyx_PyGen_FetchStopIterationValue(pvalue) \
    __Pyx_PyGen__FetchStopIterationValue($local_tstate_cname, pvalue)
#else
#define __Pyx_PyGen_FetchStopIterationValue(pvalue) \
    __Pyx_PyGen__FetchStopIterationValue(__Pyx_PyThreadState_Current, pvalue)
#endif
static int __Pyx_PyGen__FetchStopIterationValue(PyThreadState *tstate, PyObject **pvalue); /*proto*/
static CYTHON_INLINE void __Pyx_Coroutine_ResetFrameBackpointer(__Pyx_ExcInfoStruct *exc_state); /*proto*/

static char __Pyx_Coroutine_test_and_set_is_running(__pyx_CoroutineObject *gen);
static void __Pyx_Coroutine_unset_is_running(__pyx_CoroutineObject *gen);
static char __Pyx_Coroutine_get_is_running(__pyx_CoroutineObject *gen);
static PyObject *__Pyx_Coroutine_get_is_running_getter(PyObject *gen, void *closure);

#if __PYX_HAS_PY_AM_SEND == 2
static void __Pyx_SetBackportTypeAmSend(PyTypeObject *type, __Pyx_PyAsyncMethodsStruct *static_amsend_methods, __Pyx_pyiter_sendfunc am_send); /* proto */
#endif

static PyObject *__Pyx_Coroutine_fail_reduce_ex(PyObject *self, PyObject *arg); /*proto*/

//////////////////// Coroutine.proto ////////////////////

#define __Pyx_Coroutine_USED
#define __Pyx_Coroutine_CheckExact(obj) __Pyx_IS_TYPE(obj, CGLOBAL(__pyx_CoroutineType))
// __Pyx_Coroutine_Check(obj): see override for IterableCoroutine below
#define __Pyx_Coroutine_Check(obj) __Pyx_Coroutine_CheckExact(obj)
#define __Pyx_CoroutineAwait_CheckExact(obj) __Pyx_IS_TYPE(obj, CGLOBAL(__pyx_CoroutineAwaitType))

#define __Pyx_Coroutine_New(body, code, closure, name, qualname, module_name)  \
    __Pyx__Coroutine_New(CGLOBAL(__pyx_CoroutineType), body, code, closure, name, qualname, module_name)

static int __pyx_Coroutine_init(PyObject *module); /*proto*/
static PyObject *__Pyx__Coroutine_await(PyObject *coroutine); /*proto*/

typedef struct {
    PyObject_HEAD
    PyObject *coroutine;
} __pyx_CoroutineAwaitObject;

static __Pyx_PySendResult __Pyx_CoroutineAwait_Close(__pyx_CoroutineAwaitObject *self); /*proto*/


//////////////////// Generator.proto ////////////////////

#define __Pyx_Generator_USED
#define __Pyx_Generator_CheckExact(obj) __Pyx_IS_TYPE(obj, CGLOBAL(__pyx_GeneratorType))

#define __Pyx_Generator_New(body, code, closure, name, qualname, module_name)  \
    __Pyx__Coroutine_New(CGLOBAL(__pyx_GeneratorType), body, code, closure, name, qualname, module_name)

static PyObject *__Pyx_Generator_Next(PyObject *self);
static int __pyx_Generator_init(PyObject *module); /*proto*/
static CYTHON_INLINE PyObject *__Pyx_Generator_GetInlinedResult(PyObject *self); /*proto*/


//////////////////// AsyncGen ////////////////////
//@requires: AsyncGen.c::AsyncGenerator
// -> empty, only delegates to separate file


//////////////////// CoroutineBase ////////////////////
//@substitute: naming
//@requires: ReturnWithStopIteration
//@requires: Exceptions.c::PyErrFetchRestore
//@requires: Exceptions.c::PyThreadStateGet
//@requires: Exceptions.c::SwapException
//@requires: Exceptions.c::RaiseException
//@requires: Exceptions.c::SaveResetException
//@requires: ObjectHandling.c::PyObjectCallMethod1
//@requires: ObjectHandling.c::PyObjectCallNoArg
//@requires: ObjectHandling.c::PyObjectCallOneArg
//@requires: ObjectHandling.c::PyObjectFastCall
//@requires: ObjectHandling.c::PyObjectGetAttrStr
//@requires: ObjectHandling.c::PyObjectGetAttrStrNoError
//@requires: ObjectHandling.c::IterNextPlain
//@requires: CommonStructures.c::FetchCommonType
//@requires: CommonStructures.c::CommonTypesMetaclass
//@requires: ModuleSetupCode.c::IncludeStructmemberH
//@requires: ExtensionTypes.c::CallTypeTraverse
//@requires: Synchronization.c::CriticalSections

#if !CYTHON_COMPILING_IN_LIMITED_API
#include <frameobject.h>
#if PY_VERSION_HEX >= 0x030b00a6 && !defined(PYPY_VERSION)
  #ifndef Py_BUILD_CORE
    #define Py_BUILD_CORE 1
  #endif
  #include "internal/pycore_frame.h"
#endif
#endif // CYTHON_COMPILING_IN_LIMITED_API

static CYTHON_INLINE void
__Pyx_Coroutine_Undelegate(__pyx_CoroutineObject *gen) {
#if CYTHON_USE_AM_SEND
    gen->yieldfrom_am_send = NULL;
#endif
    Py_CLEAR(gen->yieldfrom);
}

//   If StopIteration exception is set, fetches its 'value'
//   attribute if any, otherwise sets pvalue to None.
//
//   Returns 0 if no exception or StopIteration is set.
//   If any other exception is set, returns -1 and leaves
//   pvalue unchanged.
static int __Pyx_PyGen__FetchStopIterationValue(PyThreadState *$local_tstate_cname, PyObject **pvalue) {
    PyObject *et, *ev, *tb;
    PyObject *value = NULL;
    CYTHON_UNUSED_VAR($local_tstate_cname);

    __Pyx_ErrFetch(&et, &ev, &tb);

    if (!et) {
        Py_XDECREF(tb);
        Py_XDECREF(ev);
        Py_INCREF(Py_None);
        *pvalue = Py_None;
        return 0;
    }

    // most common case: plain StopIteration without or with separate argument
    if (likely(et == PyExc_StopIteration)) {
        if (!ev) {
            Py_INCREF(Py_None);
            value = Py_None;
        }
        else if (likely(__Pyx_IS_TYPE(ev, (PyTypeObject*)PyExc_StopIteration))) {
            #if CYTHON_COMPILING_IN_LIMITED_API || CYTHON_COMPILING_IN_GRAAL
            value = PyObject_GetAttr(ev, PYIDENT("value"));
            if (unlikely(!value)) goto limited_api_failure;
            #else
            value = ((PyStopIterationObject *)ev)->value;
            Py_INCREF(value);
            #endif
            Py_DECREF(ev);
        }
        // PyErr_SetObject() and friends put the value directly into ev
        else if (unlikely(PyTuple_Check(ev))) {
            // if it's a tuple, it is interpreted as separate constructor arguments (surprise!)
            Py_ssize_t tuple_size = __Pyx_PyTuple_GET_SIZE(ev);
            #if !CYTHON_ASSUME_SAFE_SIZE
            if (unlikely(tuple_size < 0)) {
                Py_XDECREF(tb);
                Py_DECREF(ev);
                Py_DECREF(et);
                return -1;
            }
            #endif
            if (tuple_size >= 1) {
#if CYTHON_ASSUME_SAFE_MACROS && !CYTHON_AVOID_BORROWED_REFS
                value = PyTuple_GET_ITEM(ev, 0);
                Py_INCREF(value);
#elif CYTHON_ASSUME_SAFE_MACROS
                value = PySequence_ITEM(ev, 0);
#else
                value = PySequence_GetItem(ev, 0);
                if (!value) goto limited_api_failure;
#endif
            } else {
                Py_INCREF(Py_None);
                value = Py_None;
            }
            Py_DECREF(ev);
        }
        else if (!__Pyx_TypeCheck(ev, (PyTypeObject*)PyExc_StopIteration)) {
            // 'steal' reference to ev
            value = ev;
        }
        if (likely(value)) {
            Py_XDECREF(tb);
            Py_DECREF(et);
            *pvalue = value;
            return 0;
        }
    } else if (!__Pyx_PyErr_GivenExceptionMatches(et, PyExc_StopIteration)) {
        __Pyx_ErrRestore(et, ev, tb);
        return -1;
    }

    // otherwise: normalise and check what that gives us
    PyErr_NormalizeException(&et, &ev, &tb);
    if (unlikely(!PyObject_TypeCheck(ev, (PyTypeObject*)PyExc_StopIteration))) {
        // looks like normalisation failed - raise the new exception
        __Pyx_ErrRestore(et, ev, tb);
        return -1;
    }
    Py_XDECREF(tb);
    Py_DECREF(et);
#if CYTHON_COMPILING_IN_LIMITED_API
    value = PyObject_GetAttr(ev, PYIDENT("value"));
#else
    value = ((PyStopIterationObject *)ev)->value;
    Py_INCREF(value);
#endif
    Py_DECREF(ev);
#if CYTHON_COMPILING_IN_LIMITED_API
    if (unlikely(!value)) return -1;
#endif
    *pvalue = value;
    return 0;

#if CYTHON_COMPILING_IN_LIMITED_API || CYTHON_COMPILING_IN_GRAAL || !CYTHON_ASSUME_SAFE_MACROS
  limited_api_failure:
    Py_XDECREF(et);
    Py_XDECREF(tb);
    Py_XDECREF(ev);
    return -1;
#endif
}

static CYTHON_INLINE
__Pyx_PySendResult __Pyx_Coroutine_status_from_result(PyObject **retval) {
    if (*retval) {
        return PYGEN_NEXT;
    } else if (likely(__Pyx_PyGen__FetchStopIterationValue(__Pyx_PyThreadState_Current, retval) == 0)) {
        //assert (*retval == Py_None);
        return PYGEN_RETURN;
    } else {
        return PYGEN_ERROR;
    }
}

static CYTHON_INLINE
void __Pyx_Coroutine_ExceptionClear(__Pyx_ExcInfoStruct *exc_state) {
#if PY_VERSION_HEX >= 0x030B00a4
    Py_CLEAR(exc_state->exc_value);
#else
    PyObject *t, *v, *tb;
    t = exc_state->exc_type;
    v = exc_state->exc_value;
    tb = exc_state->exc_traceback;

    exc_state->exc_type = NULL;
    exc_state->exc_value = NULL;
    exc_state->exc_traceback = NULL;

    Py_XDECREF(t);
    Py_XDECREF(v);
    Py_XDECREF(tb);
#endif
}

#define __Pyx_Coroutine_AlreadyRunningError(gen)  (__Pyx__Coroutine_AlreadyRunningError(gen), (PyObject*)NULL)
static void __Pyx__Coroutine_AlreadyRunningError(__pyx_CoroutineObject *gen) {
    const char *msg;
    CYTHON_MAYBE_UNUSED_VAR(gen);
    if ((0)) {
    #ifdef __Pyx_Coroutine_USED
    } else if (__Pyx_Coroutine_Check((PyObject*)gen)) {
        msg = "coroutine already executing";
    #endif
    #ifdef __Pyx_AsyncGen_USED
    } else if (__Pyx_AsyncGen_CheckExact((PyObject*)gen)) {
        msg = "async generator already executing";
    #endif
    } else {
        msg = "generator already executing";
    }
    PyErr_SetString(PyExc_ValueError, msg);
}

static void __Pyx_Coroutine_AlreadyTerminatedError(PyObject *gen, PyObject *value, int closing) {
    CYTHON_MAYBE_UNUSED_VAR(gen);
    CYTHON_MAYBE_UNUSED_VAR(closing);
    #ifdef __Pyx_Coroutine_USED
    if (!closing && __Pyx_Coroutine_Check(gen)) {
        // `self` is an exhausted coroutine: raise an error,
        // except when called from gen_close(), which should
        // always be a silent method.
        PyErr_SetString(PyExc_RuntimeError, "cannot reuse already awaited coroutine");
    } else
    #endif
    if (value) {
        // `gen` is an exhausted generator:
        // only set exception if called from send().
        #ifdef __Pyx_AsyncGen_USED
        if (__Pyx_AsyncGen_CheckExact(gen))
            PyErr_SetNone(PyExc_StopAsyncIteration);
        else
        #endif
        PyErr_SetNone(PyExc_StopIteration);
    }
}

static
__Pyx_PySendResult __Pyx_Coroutine_SendEx(__pyx_CoroutineObject *self, PyObject *value, PyObject **result, int closing) {
    __Pyx_PyThreadState_declare
    PyThreadState *tstate;
    __Pyx_ExcInfoStruct *exc_state;
    PyObject *retval;

    assert(__Pyx_Coroutine_get_is_running(self));  // Callers should ensure is_running

    if (unlikely(self->resume_label == -1)) {
        __Pyx_Coroutine_AlreadyTerminatedError((PyObject*)self, value, closing);
        return PYGEN_ERROR;
    }

#if CYTHON_FAST_THREAD_STATE
    __Pyx_PyThreadState_assign
    tstate = $local_tstate_cname;
#else
    tstate = __Pyx_PyThreadState_Current;
#endif

    // Traceback/Frame rules pre-Py3.7:
    // - on entry, save external exception state in self->gi_exc_state, restore it on exit
    // - on exit, keep internally generated exceptions in self->gi_exc_state, clear everything else
    // - on entry, set "f_back" pointer of internal exception traceback to (current) outer call frame
    // - on exit, clear "f_back" of internal exception traceback
    // - do not touch external frames and tracebacks

    // Traceback/Frame rules for Py3.7+ (CYTHON_USE_EXC_INFO_STACK):
    // - on entry, push internal exception state in self->gi_exc_state on the exception stack
    // - on exit, keep internally generated exceptions in self->gi_exc_state, clear everything else
    // - on entry, set "f_back" pointer of internal exception traceback to (current) outer call frame
    // - on exit, clear "f_back" of internal exception traceback
    // - do not touch external frames and tracebacks

    // In the Limited API/PyPy
    // - on entry, save external exception state in self->gi_exc_state, restore it on exit
    // - on exit, keep internally generated exceptions in self->gi_exc_state, clear everything else
    // - don't mess with f_back because we can't work out how to do that
    // - do not touch external frames and tracebacks

    exc_state = &self->gi_exc_state;
    if (exc_state->exc_value) {
        #if CYTHON_COMPILING_IN_LIMITED_API || CYTHON_COMPILING_IN_PYPY
        // FIXME - it'd be nice to handle f_back
        #else
        // Generators always return to their most recent caller, not
        // necessarily their creator.
        PyObject *exc_tb;
        #if PY_VERSION_HEX >= 0x030B00a4 && !CYTHON_COMPILING_IN_CPYTHON
        // owned reference!
        exc_tb = PyException_GetTraceback(exc_state->exc_value);
        #elif PY_VERSION_HEX >= 0x030B00a4
        exc_tb = ((PyBaseExceptionObject*) exc_state->exc_value)->traceback;
        #else
        exc_tb = exc_state->exc_traceback;
        #endif
        if (exc_tb) {
            PyTracebackObject *tb = (PyTracebackObject *) exc_tb;
            PyFrameObject *f = tb->tb_frame;

            assert(f->f_back == NULL);
            #if PY_VERSION_HEX >= 0x030B00A1
            // PyThreadState_GetFrame returns NULL if there isn't a current frame
            // which is a valid state so no need to check
            f->f_back = PyThreadState_GetFrame(tstate);
            #else
            Py_XINCREF(tstate->frame);
            f->f_back = tstate->frame;
            #endif
            #if PY_VERSION_HEX >= 0x030B00a4 && !CYTHON_COMPILING_IN_CPYTHON
            Py_DECREF(exc_tb);
            #endif
        }
        #endif
    }

#if CYTHON_USE_EXC_INFO_STACK
    // See  https://bugs.python.org/issue25612
    exc_state->previous_item = tstate->exc_info;
    tstate->exc_info = exc_state;
#else
    if (exc_state->exc_type) {
        // We were in an except handler when we left,
        // restore the exception state which was put aside.
        __Pyx_ExceptionSwap(&exc_state->exc_type, &exc_state->exc_value, &exc_state->exc_traceback);
        // self->exc_* now holds the exception state of the caller
    } else {
        // save away the exception state of the caller
        __Pyx_Coroutine_ExceptionClear(exc_state);
        __Pyx_ExceptionSave(&exc_state->exc_type, &exc_state->exc_value, &exc_state->exc_traceback);
    }
#endif

    retval = self->body(self, tstate, value);

#if CYTHON_USE_EXC_INFO_STACK
    // See  https://bugs.python.org/issue25612
    exc_state = &self->gi_exc_state;
    tstate->exc_info = exc_state->previous_item;
    exc_state->previous_item = NULL;
    // Cut off the exception frame chain so that we can reconnect it on re-entry above.
    __Pyx_Coroutine_ResetFrameBackpointer(exc_state);
#endif

    *result = retval;

    if (self->resume_label == -1) {
        // terminated => return or error?
        return likely(retval) ? PYGEN_RETURN : PYGEN_ERROR;
    }
    return PYGEN_NEXT;
}

static CYTHON_INLINE void __Pyx_Coroutine_ResetFrameBackpointer(__Pyx_ExcInfoStruct *exc_state) {
    // Don't keep the reference to f_back any longer than necessary.  It
    // may keep a chain of frames alive or it could create a reference
    // cycle.
#if CYTHON_COMPILING_IN_PYPY || CYTHON_COMPILING_IN_LIMITED_API
    // FIXME: what to do in PyPy?
    CYTHON_UNUSED_VAR(exc_state);
#else
    PyObject *exc_tb;

    #if PY_VERSION_HEX >= 0x030B00a4
    if (!exc_state->exc_value) return;
    // owned reference!
    exc_tb = PyException_GetTraceback(exc_state->exc_value);
    #else
    exc_tb = exc_state->exc_traceback;
    #endif

    if (likely(exc_tb)) {
        PyTracebackObject *tb = (PyTracebackObject *) exc_tb;
        PyFrameObject *f = tb->tb_frame;
        Py_CLEAR(f->f_back);
        #if PY_VERSION_HEX >= 0x030B00a4
        Py_DECREF(exc_tb);
        #endif
    }
#endif
}

#define __Pyx_Coroutine_MethodReturnFromResult(gen, result, retval, iternext) \
    ((result) == PYGEN_NEXT ? (retval) : __Pyx__Coroutine_MethodReturnFromResult(gen, result, retval, iternext))

static PyObject *
__Pyx__Coroutine_MethodReturnFromResult(PyObject* gen, __Pyx_PySendResult result, PyObject *retval, int iternext) {
    CYTHON_MAYBE_UNUSED_VAR(gen);
    if (likely(result == PYGEN_RETURN)) {
        // return values pass through StopIteration or StopAsyncIteration
        int is_async = 0;
        #ifdef __Pyx_AsyncGen_USED
        is_async = __Pyx_AsyncGen_CheckExact(gen);
        #endif
        __Pyx_ReturnWithStopIteration(retval, is_async, iternext);
        Py_XDECREF(retval);
    }
    return NULL;
}

#if CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE
PyObject *__Pyx_PyGen_Send(PyGenObject *gen, PyObject *arg) {
#if PY_VERSION_HEX <= 0x030A00A1
    return _PyGen_Send(gen, arg);
#else
    PyObject *result;
    // PyIter_Send() asserts non-NULL arg
    if (PyIter_Send((PyObject*)gen, arg ? arg : Py_None, &result) == PYGEN_RETURN) {
        if (PyAsyncGen_CheckExact(gen)) {
            assert(result == Py_None);
            PyErr_SetNone(PyExc_StopAsyncIteration);
        }
        else if (result == Py_None) {
            PyErr_SetNone(PyExc_StopIteration);
        }
        else {
#if PY_VERSION_HEX < 0x030d00A1
            _PyGen_SetStopIterationValue(result);
#else
            if (!PyTuple_Check(result) && !PyExceptionInstance_Check(result)) {
                // delay instantiation if possible
                PyErr_SetObject(PyExc_StopIteration, result);
            } else {
                PyObject *exc = __Pyx_PyObject_CallOneArg(PyExc_StopIteration, result);
                if (likely(exc != NULL)) {
                    PyErr_SetObject(PyExc_StopIteration, exc);
                    Py_DECREF(exc);
                }
            }
#endif
        }
        Py_DECREF(result);
        result = NULL;
    }
    return result;
#endif
}
#endif

static CYTHON_INLINE __Pyx_PySendResult
__Pyx_Coroutine_FinishDelegation(__pyx_CoroutineObject *gen, PyObject** retval) {
    __Pyx_PySendResult result;
    PyObject *val = NULL;
    assert(__Pyx_Coroutine_get_is_running(gen));
    __Pyx_Coroutine_Undelegate(gen);
    __Pyx_PyGen__FetchStopIterationValue(__Pyx_PyThreadState_Current, &val);
    // val == NULL on failure => pass on exception
    result = __Pyx_Coroutine_SendEx(gen, val, retval, 0);
    Py_XDECREF(val);
    return result;
}

#if CYTHON_USE_AM_SEND
static __Pyx_PySendResult
__Pyx_Coroutine_SendToDelegate(__pyx_CoroutineObject *gen, __Pyx_pyiter_sendfunc gen_am_send, PyObject *value, PyObject **retval) {
    PyObject *ret = NULL;
    __Pyx_PySendResult delegate_result, result;
    assert(__Pyx_Coroutine_get_is_running(gen));
    // we assume that gen->yieldfrom cannot change as long as 'gen->is_running' is set => no safety INCREF()
    delegate_result = gen_am_send(gen->yieldfrom, value, &ret);
    if (delegate_result == PYGEN_NEXT) {
        assert (ret != NULL);
        *retval = ret;
        return PYGEN_NEXT;
    }
    //assert (delegate_result == PYGEN_ERROR ? ret == NULL : ret != NULL);
    assert (delegate_result != PYGEN_ERROR || ret == NULL);
    __Pyx_Coroutine_Undelegate(gen);
    result = __Pyx_Coroutine_SendEx(gen, ret, retval, 0);
    Py_XDECREF(ret);
    return result;
}
#endif

static PyObject *__Pyx_Coroutine_Send(PyObject *self, PyObject *value) {
    PyObject *retval = NULL;
    __Pyx_PySendResult result = __Pyx_Coroutine_AmSend(self, value, &retval);
    return __Pyx_Coroutine_MethodReturnFromResult(self, result, retval, 0);
}

static __Pyx_PySendResult
__Pyx_Coroutine_AmSend(PyObject *self, PyObject *value, PyObject **retval) {
    __Pyx_PySendResult result;
    __pyx_CoroutineObject *gen = (__pyx_CoroutineObject*) self;
    if (unlikely(__Pyx_Coroutine_test_and_set_is_running(gen))) {
        *retval = __Pyx_Coroutine_AlreadyRunningError(gen);
        return PYGEN_ERROR;
    }
    #if CYTHON_USE_AM_SEND
    if (gen->yieldfrom_am_send) {
        result = __Pyx_Coroutine_SendToDelegate(gen, gen->yieldfrom_am_send, value, retval);
    } else
    #endif
    if (gen->yieldfrom) {
        PyObject *yf = gen->yieldfrom;
        PyObject *ret;
      #if !CYTHON_USE_AM_SEND
        // Py3.10 puts "am_send" into "gen->yieldfrom_am_send" instead of using these special cases.
        // See "__Pyx_Coroutine_Set_Owned_Yield_From()" above.
        #ifdef __Pyx_Generator_USED
        if (__Pyx_Generator_CheckExact(yf)) {
            ret = __Pyx_Coroutine_Send(yf, value);
        } else
        #endif
        #ifdef __Pyx_Coroutine_USED
        if (__Pyx_Coroutine_Check(yf)) {
            ret = __Pyx_Coroutine_Send(yf, value);
        } else
        #endif
        #ifdef __Pyx_AsyncGen_USED
        if (__pyx_PyAsyncGenASend_CheckExact(yf)) {
            ret = __Pyx_async_gen_asend_send(yf, value);
        } else
        #endif
        #if CYTHON_COMPILING_IN_CPYTHON
        if (PyGen_CheckExact(yf)) {
            ret = __Pyx_PyGen_Send((PyGenObject*)yf, value == Py_None ? NULL : value);
        } else
        if (PyCoro_CheckExact(yf)) {
            ret = __Pyx_PyGen_Send((PyGenObject*)yf, value == Py_None ? NULL : value);
        } else
        #endif
      #endif
        {
            #if !CYTHON_COMPILING_IN_LIMITED_API || __PYX_LIMITED_VERSION_HEX >= 0x03080000
            // PyIter_Check() is needed here but broken in the Py3.7 Limited API.
            if (value == Py_None && PyIter_Check(yf))
                ret = __Pyx_PyIter_Next_Plain(yf);
            else
            #endif
                ret = __Pyx_PyObject_CallMethod1(yf, PYIDENT("send"), value);
        }
        if (likely(ret)) {
            __Pyx_Coroutine_unset_is_running(gen);
            *retval = ret;
            return PYGEN_NEXT;
        }
        result = __Pyx_Coroutine_FinishDelegation(gen, retval);
    } else {
        result = __Pyx_Coroutine_SendEx(gen, value, retval, 0);
    }
    __Pyx_Coroutine_unset_is_running(gen);
    return result;
}

//   This helper function is used by gen_close and gen_throw to
//   close a subiterator being delegated to by yield-from.
static int __Pyx_Coroutine_CloseIter(__pyx_CoroutineObject *gen, PyObject *yf) {
    __Pyx_PySendResult result;
    PyObject *retval = NULL;
    CYTHON_UNUSED_VAR(gen);

    assert(__Pyx_Coroutine_get_is_running(gen));

    #ifdef __Pyx_Generator_USED
    if (__Pyx_Generator_CheckExact(yf)) {
        result = __Pyx_Coroutine_Close(yf, &retval);
    } else
    #endif
    #ifdef __Pyx_Coroutine_USED
    if (__Pyx_Coroutine_Check(yf)) {
        result = __Pyx_Coroutine_Close(yf, &retval);
    } else
    if (__Pyx_CoroutineAwait_CheckExact(yf)) {
        result = __Pyx_CoroutineAwait_Close((__pyx_CoroutineAwaitObject*)yf);
    } else
    #endif
    #ifdef __Pyx_AsyncGen_USED
    if (__pyx_PyAsyncGenASend_CheckExact(yf)) {
        retval = __Pyx_async_gen_asend_close(yf, NULL);
        // cannot fail
        result = PYGEN_RETURN;
    } else
    if (__pyx_PyAsyncGenAThrow_CheckExact(yf)) {
        retval = __Pyx_async_gen_athrow_close(yf, NULL);
        // cannot fail
        result = PYGEN_RETURN;
    } else
    #endif
    {
        PyObject *meth;
        result = PYGEN_RETURN;
        meth = __Pyx_PyObject_GetAttrStrNoError(yf, PYIDENT("close"));
        if (unlikely(!meth)) {
            if (unlikely(PyErr_Occurred())) {
                PyErr_WriteUnraisable(yf);
            }
        } else {
            retval = __Pyx_PyObject_CallNoArg(meth);
            Py_DECREF(meth);
            if (unlikely(!retval)) {
                result = PYGEN_ERROR;
            }
        }
    }
    Py_XDECREF(retval);
    return result == PYGEN_ERROR ? -1 : 0;
}

static PyObject *__Pyx_Generator_Next(PyObject *self) {
    __Pyx_PySendResult result;
    PyObject *retval = NULL;
    __pyx_CoroutineObject *gen = (__pyx_CoroutineObject*) self;
    if (unlikely(__Pyx_Coroutine_test_and_set_is_running(gen))) {
        return __Pyx_Coroutine_AlreadyRunningError(gen);
    }
    #if CYTHON_USE_AM_SEND
    if (gen->yieldfrom_am_send) {
        result = __Pyx_Coroutine_SendToDelegate(gen, gen->yieldfrom_am_send, Py_None, &retval);
    } else
    #endif
    if (gen->yieldfrom) {
        PyObject *yf = gen->yieldfrom;
        PyObject *ret;
        // YieldFrom code ensures that yf is an iterator
        #ifdef __Pyx_Generator_USED
        if (__Pyx_Generator_CheckExact(yf)) {
            ret = __Pyx_Generator_Next(yf);
        } else
        #endif
        #ifdef __Pyx_Coroutine_USED
        // See https://github.com/cython/cython/issues/1999
        if (__Pyx_Coroutine_CheckExact(yf)) {
            ret = __Pyx_Coroutine_Send(yf, Py_None);
        } else
        #endif
        #if CYTHON_COMPILING_IN_CPYTHON && (PY_VERSION_HEX < 0x030A00A3 || !CYTHON_USE_AM_SEND)
        // _PyGen_Send() is not needed in 3.10+ due to "am_send"
        if (PyGen_CheckExact(yf)) {
            ret = __Pyx_PyGen_Send((PyGenObject*)yf, NULL);
        } else
        #endif
            ret = __Pyx_PyIter_Next_Plain(yf);
        if (likely(ret)) {
            __Pyx_Coroutine_unset_is_running(gen);
            return ret;
        }
        result = __Pyx_Coroutine_FinishDelegation(gen, &retval);
    } else {
        result = __Pyx_Coroutine_SendEx(gen, Py_None, &retval, 0);
    }
    __Pyx_Coroutine_unset_is_running(gen);

    return __Pyx_Coroutine_MethodReturnFromResult(self, result, retval, 1);
}

static PyObject *__Pyx_Coroutine_Close_Method(PyObject *self, PyObject *arg) {
    PyObject *retval = NULL;
    __Pyx_PySendResult result;
    CYTHON_UNUSED_VAR(arg);
    result = __Pyx_Coroutine_Close(self, &retval);
    if (unlikely(result == PYGEN_ERROR))
        return NULL;
    Py_XDECREF(retval);
    Py_RETURN_NONE;
}

static __Pyx_PySendResult
__Pyx_Coroutine_Close(PyObject *self, PyObject **retval) {
    __pyx_CoroutineObject *gen = (__pyx_CoroutineObject *) self;
    __Pyx_PySendResult result;
    PyObject *yf;
    int err = 0;

    if (unlikely(__Pyx_Coroutine_test_and_set_is_running(gen))) {
        *retval = __Pyx_Coroutine_AlreadyRunningError(gen);
        return PYGEN_ERROR;
    }
    yf = gen->yieldfrom;

    if (yf) {
        Py_INCREF(yf);
        err = __Pyx_Coroutine_CloseIter(gen, yf);
        __Pyx_Coroutine_Undelegate(gen);
        Py_DECREF(yf);
    }
    if (err == 0)
        PyErr_SetNone(PyExc_GeneratorExit);
    result = __Pyx_Coroutine_SendEx(gen, NULL, retval, 1);
    if (result == PYGEN_ERROR) {
        // WARNING: *retval == NULL !
        __Pyx_PyThreadState_declare
        __Pyx_PyThreadState_assign
        __Pyx_Coroutine_unset_is_running(gen);
        if (!__Pyx_PyErr_Occurred()) {
            return PYGEN_RETURN;
        } else if (likely(__Pyx_PyErr_ExceptionMatches2(PyExc_GeneratorExit, PyExc_StopIteration))) {
            __Pyx_PyErr_Clear();
            return PYGEN_RETURN;
        }
        return PYGEN_ERROR;
    } else if (likely(result == PYGEN_RETURN && *retval == Py_None)) {
        __Pyx_Coroutine_unset_is_running(gen);
        return PYGEN_RETURN;
    } else {
        const char *msg;
        Py_DECREF(*retval);
        *retval = NULL;
        if ((0)) {
        #ifdef __Pyx_Coroutine_USED
        } else if (__Pyx_Coroutine_Check(self)) {
            msg = "coroutine ignored GeneratorExit";
        #endif
        #ifdef __Pyx_AsyncGen_USED
        } else if (__Pyx_AsyncGen_CheckExact(self)) {
            msg = "async generator ignored GeneratorExit";
        #endif
        } else {
            msg = "generator ignored GeneratorExit";
        }
        PyErr_SetString(PyExc_RuntimeError, msg);
        __Pyx_Coroutine_unset_is_running(gen);
        return PYGEN_ERROR;
    }
}

static PyObject *__Pyx__Coroutine_Throw(PyObject *self, PyObject *typ, PyObject *val, PyObject *tb,
                                        PyObject *args, int close_on_genexit) {
    __pyx_CoroutineObject *gen = (__pyx_CoroutineObject *) self;
    PyObject *yf;

    if (unlikely(__Pyx_Coroutine_test_and_set_is_running(gen)))
        return __Pyx_Coroutine_AlreadyRunningError(gen);

    yf = gen->yieldfrom;
    if (yf) {
        __Pyx_PySendResult result;
        PyObject *ret;
        Py_INCREF(yf);
        if (__Pyx_PyErr_GivenExceptionMatches(typ, PyExc_GeneratorExit) && close_on_genexit) {
            // Asynchronous generators *should not* be closed right away.
            // We have to allow some awaits to work it through, hence the
            // `close_on_genexit` parameter here.
            int err = __Pyx_Coroutine_CloseIter(gen, yf);
            Py_DECREF(yf);
            __Pyx_Coroutine_Undelegate(gen);
            if (err < 0)
                goto propagate_exception;
            goto throw_here;
        }
        if (0
        #ifdef __Pyx_Generator_USED
            || __Pyx_Generator_CheckExact(yf)
        #endif
        #ifdef __Pyx_Coroutine_USED
            || __Pyx_Coroutine_Check(yf)
        #endif
            ) {
            ret = __Pyx__Coroutine_Throw(yf, typ, val, tb, args, close_on_genexit);
        #ifdef __Pyx_Coroutine_USED
        } else if (__Pyx_CoroutineAwait_CheckExact(yf)) {
            ret = __Pyx__Coroutine_Throw(((__pyx_CoroutineAwaitObject*)yf)->coroutine, typ, val, tb, args, close_on_genexit);
        #endif
        } else {
            PyObject *meth = __Pyx_PyObject_GetAttrStrNoError(yf, PYIDENT("throw"));
            if (unlikely(!meth)) {
                Py_DECREF(yf);
                if (unlikely(PyErr_Occurred())) {
                    __Pyx_Coroutine_unset_is_running(gen);
                    return NULL;
                }
                __Pyx_Coroutine_Undelegate(gen);
                goto throw_here;
            }
            if (likely(args)) {
                ret = __Pyx_PyObject_Call(meth, args, NULL);
            } else {
                // "tb" or even "val" might be NULL.
                PyObject *cargs[4] = {NULL, typ, val, tb};
                size_t nargs = 1U + (val != NULL) + (tb != NULL);
                ret = __Pyx_PyObject_FastCall(meth, cargs+1, nargs | __Pyx_PY_VECTORCALL_ARGUMENTS_OFFSET);
            }
            Py_DECREF(meth);
        }
        Py_DECREF(yf);
        if (ret) {
            __Pyx_Coroutine_unset_is_running(gen);
            return ret;
        }

        result = __Pyx_Coroutine_FinishDelegation(gen, &ret);
        __Pyx_Coroutine_unset_is_running(gen);
        return __Pyx_Coroutine_MethodReturnFromResult(self, result, ret, 0);
    }
throw_here:
    __Pyx_Raise(typ, val, tb, NULL);

propagate_exception:
    {
        PyObject *retval = NULL;
        __Pyx_PySendResult result = __Pyx_Coroutine_SendEx(gen, NULL, &retval, 0);
        __Pyx_Coroutine_unset_is_running(gen);
        return __Pyx_Coroutine_MethodReturnFromResult(self, result, retval, 0);
    }
}

static PyObject *__Pyx_Coroutine_Throw(PyObject *self, PyObject *args) {
    PyObject *typ;
    PyObject *val = NULL;
    PyObject *tb = NULL;

    if (unlikely(!PyArg_UnpackTuple(args, "throw", 1, 3, &typ, &val, &tb)))
        return NULL;

    return __Pyx__Coroutine_Throw(self, typ, val, tb, args, 1);
}

static CYTHON_INLINE int __Pyx_Coroutine_traverse_excstate(__Pyx_ExcInfoStruct *exc_state, visitproc visit, void *arg) {
#if PY_VERSION_HEX >= 0x030B00a4
    Py_VISIT(exc_state->exc_value);
#else
    Py_VISIT(exc_state->exc_type);
    Py_VISIT(exc_state->exc_value);
    Py_VISIT(exc_state->exc_traceback);
#endif
    return 0;
}

static int __Pyx_Coroutine_traverse(__pyx_CoroutineObject *gen, visitproc visit, void *arg) {
    {
        int e = __Pyx_call_type_traverse((PyObject*)gen, 1, visit, arg);
        if (e) return e;
    }
    Py_VISIT(gen->closure);
    Py_VISIT(gen->classobj);
    Py_VISIT(gen->yieldfrom);
    return __Pyx_Coroutine_traverse_excstate(&gen->gi_exc_state, visit, arg);
}

static int __Pyx_Coroutine_clear(PyObject *self) {
    __pyx_CoroutineObject *gen = (__pyx_CoroutineObject *) self;

    Py_CLEAR(gen->closure);
    Py_CLEAR(gen->classobj);
    __Pyx_Coroutine_Undelegate(gen);
    __Pyx_Coroutine_ExceptionClear(&gen->gi_exc_state);
#ifdef __Pyx_AsyncGen_USED
    if (__Pyx_AsyncGen_CheckExact(self)) {
        Py_CLEAR(((__pyx_PyAsyncGenObject*)gen)->ag_finalizer);
    }
#endif
    Py_CLEAR(gen->gi_code);
    Py_CLEAR(gen->gi_frame);
    Py_CLEAR(gen->gi_name);
    Py_CLEAR(gen->gi_qualname);
    Py_CLEAR(gen->gi_modulename);
    return 0;
}

static void __Pyx_Coroutine_dealloc(PyObject *self) {
    __pyx_CoroutineObject *gen = (__pyx_CoroutineObject *) self;

    PyObject_GC_UnTrack(gen);

    #if PY_VERSION_HEX < 0x030C0000 || CYTHON_COMPILING_IN_LIMITED_API
    if (gen->gi_weakreflist != NULL)
    #endif
        PyObject_ClearWeakRefs(self);

    if (gen->resume_label >= 0) {
        // Generator is paused or unstarted, so we need to close
        PyObject_GC_Track(self);
#if CYTHON_USE_TP_FINALIZE
        if (unlikely(PyObject_CallFinalizerFromDealloc(self)))
#else
        {
            destructor del = __Pyx_PyObject_GetSlot(gen, tp_del, destructor);
            if (del) del(self);
        }
        if (unlikely(Py_REFCNT(self) > 0))
#endif
        {
            // resurrected.  :(
            return;
        }
        PyObject_GC_UnTrack(self);
    }

#ifdef __Pyx_AsyncGen_USED
    if (__Pyx_AsyncGen_CheckExact(self)) {
        /* We have to handle this case for asynchronous generators
           right here, because this code has to be between UNTRACK
           and GC_Del. */
        Py_CLEAR(((__pyx_PyAsyncGenObject*)self)->ag_finalizer);
    }
#endif
    __Pyx_Coroutine_clear(self);
    __Pyx_PyHeapTypeObject_GC_Del(gen);
}

#if CYTHON_USE_TP_FINALIZE
static void __Pyx_Coroutine_del(PyObject *self) {
    PyObject *error_type, *error_value, *error_traceback;
    __pyx_CoroutineObject *gen = (__pyx_CoroutineObject *) self;
    __Pyx_PyThreadState_declare

    if (gen->resume_label < 0) {
        // already terminated => nothing to clean up
        return;
    }

    __Pyx_PyThreadState_assign

    // Save the current exception, if any.
    __Pyx_ErrFetch(&error_type, &error_value, &error_traceback);

#ifdef __Pyx_AsyncGen_USED
    if (__Pyx_AsyncGen_CheckExact(self)) {
        __pyx_PyAsyncGenObject *agen = (__pyx_PyAsyncGenObject*)self;
        PyObject *finalizer = agen->ag_finalizer;
        if (finalizer && !agen->ag_closed) {
            PyObject *res = __Pyx_PyObject_CallOneArg(finalizer, self);
            if (unlikely(!res)) {
                PyErr_WriteUnraisable(self);
            } else {
                Py_DECREF(res);
            }
            // Restore the saved exception.
            __Pyx_ErrRestore(error_type, error_value, error_traceback);
            return;
        }
    }
#endif

    if (unlikely(gen->resume_label == 0 && !error_value)) {
#ifdef __Pyx_Coroutine_USED
#ifdef __Pyx_Generator_USED
    // only warn about (async) coroutines
    if (!__Pyx_Generator_CheckExact(self))
#endif
        {
        // untrack dead object as we are executing Python code (which might trigger GC)
        PyObject_GC_UnTrack(self);
        if (unlikely(PyErr_WarnFormat(PyExc_RuntimeWarning, 1, "coroutine '%.50S' was never awaited", gen->gi_qualname) < 0))
            PyErr_WriteUnraisable(self);
        PyObject_GC_Track(self);
        }
#endif /*__Pyx_Coroutine_USED*/
    } else {
        PyObject *retval = NULL;
        __Pyx_PySendResult result = __Pyx_Coroutine_Close(self, &retval);
        if (result == PYGEN_ERROR) {
            PyErr_WriteUnraisable(self);
        } else {
            Py_XDECREF(retval);
        }
    }

    // Restore the saved exception.
    __Pyx_ErrRestore(error_type, error_value, error_traceback);
}
#endif

static PyObject *
__Pyx_Coroutine_get_name(__pyx_CoroutineObject *self, void *context)
{
    PyObject *name = self->gi_name;
    CYTHON_UNUSED_VAR(context);
    // avoid NULL pointer dereference during garbage collection
    if (unlikely(!name)) name = Py_None;
    Py_INCREF(name);
    return name;
}

static int
__Pyx_Coroutine_set_name(__pyx_CoroutineObject *self, PyObject *value, void *context)
{
    CYTHON_UNUSED_VAR(context);
    if (unlikely(value == NULL || !PyUnicode_Check(value))) {
        PyErr_SetString(PyExc_TypeError,
                        "__name__ must be set to a string object");
        return -1;
    }
    Py_INCREF(value);
    __Pyx_Py_XDECREF_SET(self->gi_name, value);
    return 0;
}

static PyObject *
__Pyx_Coroutine_get_qualname(__pyx_CoroutineObject *self, void *context)
{
    PyObject *name = self->gi_qualname;
    CYTHON_UNUSED_VAR(context);
    // avoid NULL pointer dereference during garbage collection
    if (unlikely(!name)) name = Py_None;
    Py_INCREF(name);
    return name;
}

static int
__Pyx_Coroutine_set_qualname(__pyx_CoroutineObject *self, PyObject *value, void *context)
{
    CYTHON_UNUSED_VAR(context);
    if (unlikely(value == NULL || !PyUnicode_Check(value))) {
        PyErr_SetString(PyExc_TypeError,
                        "__qualname__ must be set to a string object");
        return -1;
    }
    Py_INCREF(value);
    __Pyx_Py_XDECREF_SET(self->gi_qualname, value);
    return 0;
}

#if !CYTHON_COMPILING_IN_LIMITED_API
static PyObject *
__Pyx__Coroutine_get_frame_locked(__pyx_CoroutineObject *self) {
    PyObject *frame;
    frame = self->gi_frame;
    if (!frame) {
        if (unlikely(!self->gi_code)) {
            // Avoid doing something stupid, e.g. during garbage collection.
            Py_RETURN_NONE;
        }
        // The coroutine doesn't know what module it's in so can't fill in proper globals.
        // In principle this could be solved by storing a reference to the module in the coroutine,
        // but in practice it probably isn't worth the extra memory.
        // Therefore, just supply a blank dict.
        PyObject *globals = PyDict_New();
        if (unlikely(!globals)) return NULL;
        frame = (PyObject *) PyFrame_New(
            PyThreadState_Get(),            /*PyThreadState *tstate,*/
            (PyCodeObject*) self->gi_code,  /*PyCodeObject *code,*/
            globals,                        /*PyObject *globals,*/
            0                               /*PyObject *locals*/
        );
        Py_DECREF(globals);
        if (unlikely(!frame))
            return NULL;
        // handle potential race initializing the frame
        if (unlikely(self->gi_frame)) {
            Py_DECREF(frame);
            frame = self->gi_frame;
        } else {
            // keep the frame cached once it's created
            self->gi_frame = frame;
        }
    }
    Py_INCREF(frame);
    return frame;
}
#endif


static PyObject *
__Pyx__Coroutine_get_frame(__pyx_CoroutineObject *self)
{
#if !CYTHON_COMPILING_IN_LIMITED_API
    PyObject *frame;
    __Pyx_BEGIN_CRITICAL_SECTION((PyObject*)self);
    frame = __Pyx__Coroutine_get_frame_locked(self);
    __Pyx_END_CRITICAL_SECTION();
    return frame;
#else
    // In the limited API there probably isn't much we can usefully do to get a frame
    CYTHON_UNUSED_VAR(self);
    Py_RETURN_NONE;
#endif
}

static PyObject *
__Pyx_Coroutine_get_frame(__pyx_CoroutineObject *self, void *context) {
    CYTHON_UNUSED_VAR(context);
    PyObject *frame = self->gi_frame;
    if (frame)
        return __Pyx_NewRef(frame);
    return __Pyx__Coroutine_get_frame(self);
}

static __pyx_CoroutineObject *__Pyx__Coroutine_New(
            PyTypeObject* type, __pyx_coroutine_body_t body, PyObject *code, PyObject *closure,
            PyObject *name, PyObject *qualname, PyObject *module_name) {
    __pyx_CoroutineObject *gen = PyObject_GC_New(__pyx_CoroutineObject, type);
    if (unlikely(!gen))
        return NULL;
    return __Pyx__Coroutine_NewInit(gen, body, code, closure, name, qualname, module_name);
}

static __pyx_CoroutineObject *__Pyx__Coroutine_NewInit(
            __pyx_CoroutineObject *gen, __pyx_coroutine_body_t body, PyObject *code, PyObject *closure,
            PyObject *name, PyObject *qualname, PyObject *module_name) {
    gen->body = body;
    gen->closure = closure;
    Py_XINCREF(closure);
    gen->is_running = 0;
    gen->resume_label = 0;
    gen->classobj = NULL;
    gen->yieldfrom = NULL;
    gen->yieldfrom_am_send = NULL;
    #if PY_VERSION_HEX >= 0x030B00a4 && !CYTHON_COMPILING_IN_LIMITED_API
    gen->gi_exc_state.exc_value = NULL;
    #else
    gen->gi_exc_state.exc_type = NULL;
    gen->gi_exc_state.exc_value = NULL;
    gen->gi_exc_state.exc_traceback = NULL;
    #endif
#if CYTHON_USE_EXC_INFO_STACK
    gen->gi_exc_state.previous_item = NULL;
#endif
#if PY_VERSION_HEX < 0x030C0000 || CYTHON_COMPILING_IN_LIMITED_API
    gen->gi_weakreflist = NULL;
#endif
    Py_XINCREF(qualname);
    gen->gi_qualname = qualname;
    Py_XINCREF(name);
    gen->gi_name = name;
    Py_XINCREF(module_name);
    gen->gi_modulename = module_name;
    Py_XINCREF(code);
    gen->gi_code = code;
    gen->gi_frame = NULL;
#if CYTHON_USE_SYS_MONITORING && (CYTHON_PROFILE || CYTHON_TRACE)
    memset(gen->$monitoring_states_cname, 0, sizeof(gen->$monitoring_states_cname));
    gen->$monitoring_version_cname = 0;
#endif

    PyObject_GC_Track(gen);
    return gen;
}


static char __Pyx_Coroutine_test_and_set_is_running(__pyx_CoroutineObject *gen) {
    // Working on the basis that:
    // 1. is_running is the main thing needed to keep generators thread-safe. They can't be
    //    run simultaneously from multiple threads so as long as we track this correctly
    //    then all is good.
    // 2. They aren't *ultra* high performance, and so critical sections (locking) is appropriate
    //    instead of atomics.
    char result;
    __Pyx_BEGIN_CRITICAL_SECTION((PyObject*)gen);
    result = gen->is_running;
    gen->is_running = 1;
    __Pyx_END_CRITICAL_SECTION();
    return result;
}
static void __Pyx_Coroutine_unset_is_running(__pyx_CoroutineObject *gen) {
    __Pyx_BEGIN_CRITICAL_SECTION((PyObject*)gen);
    assert(gen->is_running);
    gen->is_running = 0;
    __Pyx_END_CRITICAL_SECTION();
}
static char __Pyx_Coroutine_get_is_running(__pyx_CoroutineObject *gen) {
    char result;
    __Pyx_BEGIN_CRITICAL_SECTION((PyObject*)gen);
    result = gen->is_running;
    __Pyx_END_CRITICAL_SECTION();
    return result;
}
static PyObject *__Pyx_Coroutine_get_is_running_getter(PyObject *gen, void *closure) {
    CYTHON_UNUSED_VAR(closure);
    char result = __Pyx_Coroutine_get_is_running((__pyx_CoroutineObject*)gen);
    if (result) Py_RETURN_TRUE;
    else Py_RETURN_FALSE;
}

#if __PYX_HAS_PY_AM_SEND == 2
static void __Pyx_SetBackportTypeAmSend(PyTypeObject *type, __Pyx_PyAsyncMethodsStruct *static_amsend_methods, __Pyx_pyiter_sendfunc am_send) {
    Py_ssize_t ptr_offset = (char*)(type->tp_as_async) - (char*)type;
    if (ptr_offset < 0 || ptr_offset > type->tp_basicsize) {
        // The pointer isn't to somewhere within the type. This must be a cached type that's already been updated.
        return;
    }

    // Copy the standard Python bits of it.
    memcpy((void*)static_amsend_methods, (void*)(type->tp_as_async), sizeof(*type->tp_as_async));
    static_amsend_methods->am_send = am_send;

    // Replace.
    type->tp_as_async = __Pyx_SlotTpAsAsync(static_amsend_methods);
}
#endif


// In earlier versions of Python an object with no __dict__ and not __slots__ is assumed
// to be pickleable by default. Coroutine-wrappers have significant state so shouldn't be.
// Therefore provide a default implementation.
// Something similar applies to heaptypes (i.e. with type_specs) with protocols 0 and 1
// even in more recent versions.
// We are applying this to all Python versions (hence the commented out version guard)
// to make the behaviour explicit.
// #if CYTHON_USE_TYPE_SPECS
static PyObject *__Pyx_Coroutine_fail_reduce_ex(PyObject *self, PyObject *arg) {
    CYTHON_UNUSED_VAR(arg);

    __Pyx_TypeName self_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE((PyObject*)self));
    PyErr_Format(PyExc_TypeError, "cannot pickle '" __Pyx_FMT_TYPENAME "' object",
                         self_type_name);
    __Pyx_DECREF_TypeName(self_type_name);
    return NULL;
}
// #endif

//////////////////// Coroutine.module_state_decls ////////////////////

PyTypeObject *__pyx_CoroutineType;
PyTypeObject *__pyx_CoroutineAwaitType;

//////////////////// Coroutine.module_state_traverse ///////////////////

Py_VISIT(traverse_module_state->__pyx_CoroutineType);
Py_VISIT(traverse_module_state->__pyx_CoroutineAwaitType);

//////////////////// Coroutine.module_state_clear ///////////////////

Py_CLEAR(clear_module_state->__pyx_CoroutineType);
Py_CLEAR(clear_module_state->__pyx_CoroutineAwaitType);

//////////////////// Coroutine.init //////////////////
//@substitute: naming

if (likely(__pyx_Coroutine_init($module_cname) == 0)); else

//////////////////// Coroutine ////////////////////
//@requires: CoroutineBase
//@requires: ExtensionTypes.c::CallTypeTraverse
//@substitute: naming

static void __Pyx_CoroutineAwait_dealloc(PyObject *self) {
    PyObject_GC_UnTrack(self);
    Py_CLEAR(((__pyx_CoroutineAwaitObject*)self)->coroutine);
    __Pyx_PyHeapTypeObject_GC_Del(self);
}

static int __Pyx_CoroutineAwait_traverse(__pyx_CoroutineAwaitObject *self, visitproc visit, void *arg) {
    {
        int e = __Pyx_call_type_traverse((PyObject*)self, 1, visit, arg);
        if (e) return e;
    }
    Py_VISIT(self->coroutine);
    return 0;
}

static int __Pyx_CoroutineAwait_clear(__pyx_CoroutineAwaitObject *self) {
    Py_CLEAR(self->coroutine);
    return 0;
}

static PyObject *__Pyx_CoroutineAwait_Next(__pyx_CoroutineAwaitObject *self) {
    return __Pyx_Generator_Next(self->coroutine);
}

static PyObject *__Pyx_CoroutineAwait_Send(__pyx_CoroutineAwaitObject *self, PyObject *value) {
    return __Pyx_Coroutine_Send(self->coroutine, value);
}

#if __PYX_HAS_PY_AM_SEND
static __Pyx_PySendResult __Pyx_CoroutineAwait_AmSend(PyObject *self, PyObject *value, PyObject **retval) {
    return __Pyx_Coroutine_AmSend(((__pyx_CoroutineAwaitObject*)self)->coroutine, value, retval);
}
#endif

static PyObject *__Pyx_CoroutineAwait_Throw(__pyx_CoroutineAwaitObject *self, PyObject *args) {
    return __Pyx_Coroutine_Throw(self->coroutine, args);
}

static __Pyx_PySendResult __Pyx_CoroutineAwait_Close(__pyx_CoroutineAwaitObject *self) {
    PyObject *retval = NULL;
    __Pyx_PySendResult result = __Pyx_Coroutine_Close(self->coroutine, &retval);
    Py_XDECREF(retval);
    return result;
}

static PyObject *__Pyx_CoroutineAwait_Close_Method(__pyx_CoroutineAwaitObject *self, PyObject *arg) {
    PyObject *retval = NULL;
    __Pyx_PySendResult result;
    CYTHON_UNUSED_VAR(arg);
    result = __Pyx_Coroutine_Close(self->coroutine, &retval);
    if (unlikely(result == PYGEN_ERROR))
        return NULL;
    Py_XDECREF(retval);
    Py_RETURN_NONE;
}

static PyObject *__Pyx_CoroutineAwait_self(PyObject *self) {
    Py_INCREF(self);
    return self;
}

#if !CYTHON_COMPILING_IN_PYPY
static PyObject *__Pyx_CoroutineAwait_no_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
    CYTHON_UNUSED_VAR(type);
    CYTHON_UNUSED_VAR(args);
    CYTHON_UNUSED_VAR(kwargs);
    PyErr_SetString(PyExc_TypeError, "cannot instantiate type, use 'await coroutine' instead");
    return NULL;
}
#endif

static PyMethodDef __pyx_CoroutineAwait_methods[] = {
    {"send", (PyCFunction) __Pyx_CoroutineAwait_Send, METH_O,
     PyDoc_STR("send(arg) -> send 'arg' into coroutine,\nreturn next yielded value or raise StopIteration.")},
    {"throw", (PyCFunction) __Pyx_CoroutineAwait_Throw, METH_VARARGS,
     PyDoc_STR("throw(typ[,val[,tb]]) -> raise exception in coroutine,\nreturn next yielded value or raise StopIteration.")},
    {"close", (PyCFunction) __Pyx_CoroutineAwait_Close_Method, METH_NOARGS, PyDoc_STR("close() -> raise GeneratorExit inside coroutine.")},
// only needed with type-specs, but included in all versions for clarity
// #if CYTHON_USE_TYPE_SPECS
    {"__reduce_ex__", (PyCFunction) __Pyx_Coroutine_fail_reduce_ex, METH_O, 0},
    {"__reduce__", (PyCFunction) __Pyx_Coroutine_fail_reduce_ex, METH_NOARGS, 0},
// #endif
    {0, 0, 0, 0}
};

static PyType_Slot __pyx_CoroutineAwaitType_slots[] = {
    {Py_tp_dealloc, (void *)__Pyx_CoroutineAwait_dealloc},
    {Py_tp_traverse, (void *)__Pyx_CoroutineAwait_traverse},
    {Py_tp_clear, (void *)__Pyx_CoroutineAwait_clear},
#if !CYTHON_COMPILING_IN_PYPY
    {Py_tp_new, (void *)__Pyx_CoroutineAwait_no_new},
#endif
    {Py_tp_methods, (void *)__pyx_CoroutineAwait_methods},
    {Py_tp_iter, (void *)__Pyx_CoroutineAwait_self},
    {Py_tp_iternext, (void *)__Pyx_CoroutineAwait_Next},
#if __PYX_HAS_PY_AM_SEND == 1
    {Py_am_send, (void *)__Pyx_CoroutineAwait_AmSend},
#endif
    {0, 0},
};

static PyType_Spec __pyx_CoroutineAwaitType_spec = {
    __PYX_TYPE_MODULE_PREFIX "coroutine_wrapper",
    sizeof(__pyx_CoroutineAwaitObject),
    0,
    Py_TPFLAGS_IMMUTABLETYPE | Py_TPFLAGS_DISALLOW_INSTANTIATION |
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | __Pyx_TPFLAGS_HAVE_AM_SEND, /*tp_flags*/
    __pyx_CoroutineAwaitType_slots
};

#if __PYX_HAS_PY_AM_SEND == 2 // backport
static __Pyx_PyAsyncMethodsStruct __pyx_CoroutineAwait_as_async;
#endif

static CYTHON_INLINE PyObject *__Pyx__Coroutine_await(PyObject *coroutine) {
    __pyx_CoroutineAwaitObject *await = PyObject_GC_New(__pyx_CoroutineAwaitObject, CGLOBAL(__pyx_CoroutineAwaitType));
    if (unlikely(!await)) return NULL;
    Py_INCREF(coroutine);
    await->coroutine = coroutine;
    PyObject_GC_Track(await);
    return (PyObject*)await;
}

static PyObject *__Pyx_Coroutine_await(PyObject *coroutine) {
    if (unlikely(!coroutine || !__Pyx_Coroutine_Check(coroutine))) {
        PyErr_SetString(PyExc_TypeError, "invalid input, expected coroutine");
        return NULL;
    }
    return __Pyx__Coroutine_await(coroutine);
}

static PyMethodDef __pyx_Coroutine_methods[] = {
    {"send", (PyCFunction) __Pyx_Coroutine_Send, METH_O,
     PyDoc_STR("send(arg) -> send 'arg' into coroutine,\nreturn next iterated value or raise StopIteration.")},
    {"throw", (PyCFunction) __Pyx_Coroutine_Throw, METH_VARARGS,
     PyDoc_STR("throw(typ[,val[,tb]]) -> raise exception in coroutine,\nreturn next iterated value or raise StopIteration.")},
    {"close", (PyCFunction) __Pyx_Coroutine_Close_Method, METH_NOARGS, PyDoc_STR("close() -> raise GeneratorExit inside coroutine.")},
    {"__reduce_ex__", (PyCFunction) __Pyx_Coroutine_fail_reduce_ex, METH_O, 0},
    {"__reduce__", (PyCFunction) __Pyx_Coroutine_fail_reduce_ex, METH_NOARGS, 0},
    {0, 0, 0, 0}
};

static PyMemberDef __pyx_Coroutine_memberlist[] = {
    {"cr_await", T_OBJECT, offsetof(__pyx_CoroutineObject, yieldfrom), READONLY,
     PyDoc_STR("object being awaited, or None")},
    {"cr_code", T_OBJECT, offsetof(__pyx_CoroutineObject, gi_code), READONLY, NULL},
    {"__module__", T_OBJECT, offsetof(__pyx_CoroutineObject, gi_modulename), 0, 0},
#if PY_VERSION_HEX < 0x030C0000 || CYTHON_COMPILING_IN_LIMITED_API
    {"__weaklistoffset__", T_PYSSIZET, offsetof(__pyx_CoroutineObject, gi_weakreflist), READONLY, 0},
#endif
    {0, 0, 0, 0, 0}
};

static PyGetSetDef __pyx_Coroutine_getsets[] = {
    {"__name__", (getter)__Pyx_Coroutine_get_name, (setter)__Pyx_Coroutine_set_name,
     PyDoc_STR("name of the coroutine"), 0},
    {"__qualname__", (getter)__Pyx_Coroutine_get_qualname, (setter)__Pyx_Coroutine_set_qualname,
     PyDoc_STR("qualified name of the coroutine"), 0},
    {"cr_frame", (getter)__Pyx_Coroutine_get_frame, NULL,
     PyDoc_STR("Frame of the coroutine"), 0},
    // getter rather than member for thread safety.
    {"cr_running", __Pyx_Coroutine_get_is_running_getter, NULL, NULL, NULL},
    {0, 0, 0, 0, 0}
};

static PyType_Slot __pyx_CoroutineType_slots[] = {
    {Py_tp_dealloc, (void *)__Pyx_Coroutine_dealloc},
    {Py_am_await, (void *)&__Pyx_Coroutine_await},
    {Py_tp_traverse, (void *)__Pyx_Coroutine_traverse},
    {Py_tp_methods, (void *)__pyx_Coroutine_methods},
    {Py_tp_members, (void *)__pyx_Coroutine_memberlist},
    {Py_tp_getset, (void *)__pyx_Coroutine_getsets},
    {Py_tp_getattro, (void *) PyObject_GenericGetAttr},
#if CYTHON_USE_TP_FINALIZE
    {Py_tp_finalize, (void *)__Pyx_Coroutine_del},
#endif
#if __PYX_HAS_PY_AM_SEND == 1
    {Py_am_send, (void *)__Pyx_Coroutine_AmSend},
#endif
    {0, 0},
};

static PyType_Spec __pyx_CoroutineType_spec = {
    __PYX_TYPE_MODULE_PREFIX "coroutine",
    sizeof(__pyx_CoroutineObject),
    0,
#if PY_VERSION_HEX >= 0x030C0000 && !CYTHON_COMPILING_IN_LIMITED_API
    Py_TPFLAGS_MANAGED_WEAKREF |
#endif
    Py_TPFLAGS_IMMUTABLETYPE | Py_TPFLAGS_DISALLOW_INSTANTIATION |
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | __Pyx_TPFLAGS_HAVE_AM_SEND, /*tp_flags*/
    __pyx_CoroutineType_slots
};

#if __PYX_HAS_PY_AM_SEND == 2
static __Pyx_PyAsyncMethodsStruct __pyx_Coroutine_as_async;
#endif

static int __pyx_Coroutine_init(PyObject *module) {
    $modulestatetype_cname *mstate;
    CYTHON_MAYBE_UNUSED_VAR(module);
    // on Windows, C-API functions can't be used in slots statically
    mstate = __Pyx_PyModule_GetState(module);
    mstate->__pyx_CoroutineType = __Pyx_FetchCommonTypeFromSpec(
        mstate->__pyx_CommonTypesMetaclassType, module, &__pyx_CoroutineType_spec, NULL);
    if (unlikely(!mstate->__pyx_CoroutineType))
        return -1;
#if __PYX_HAS_PY_AM_SEND == 2
    __Pyx_SetBackportTypeAmSend(mstate->__pyx_CoroutineType, &__pyx_Coroutine_as_async, &__Pyx_Coroutine_AmSend);
#endif

#ifdef __Pyx_IterableCoroutine_USED
    if (unlikely(__pyx_IterableCoroutine_init(module) == -1))
        return -1;
#endif

    mstate->__pyx_CoroutineAwaitType = __Pyx_FetchCommonTypeFromSpec(NULL, module, &__pyx_CoroutineAwaitType_spec, NULL);
    if (unlikely(!mstate->__pyx_CoroutineAwaitType))
        return -1;
#if __PYX_HAS_PY_AM_SEND == 2
    __Pyx_SetBackportTypeAmSend(mstate->__pyx_CoroutineAwaitType, &__pyx_CoroutineAwait_as_async, &__Pyx_CoroutineAwait_AmSend);
#endif

    return 0;
}

//////////////////// IterableCoroutine.module_state_decls ////////////////////

PyTypeObject *__pyx_IterableCoroutineType;

//////////////////// IterableCoroutine.module_state_traverse ///////////////////

Py_VISIT(traverse_module_state->__pyx_IterableCoroutineType);

//////////////////// IterableCoroutine.module_state_clear ///////////////////

Py_CLEAR(clear_module_state->__pyx_IterableCoroutineType);

//////////////////// IterableCoroutine.init //////////////////
//@substitute: naming

if (likely(__pyx_IterableCoroutine_init($module_cname) == 0)); else

//////////////////// IterableCoroutine.proto ////////////////////

#define __Pyx_IterableCoroutine_USED

#undef __Pyx_Coroutine_Check
#define __Pyx_Coroutine_Check(obj) (__Pyx_Coroutine_CheckExact(obj) || __Pyx_IS_TYPE(obj, CGLOBAL(__pyx_IterableCoroutineType)))

#define __Pyx_IterableCoroutine_New(body, code, closure, name, qualname, module_name)  \
    __Pyx__Coroutine_New(CGLOBAL(__pyx_IterableCoroutineType), body, code, closure, name, qualname, module_name)

static int __pyx_IterableCoroutine_init(PyObject *module);/*proto*/


//////////////////// IterableCoroutine ////////////////////
//@requires: Coroutine
//@requires: CommonStructures.c::FetchCommonType
//@requires: CommonStructures.c::CommonTypesMetaclass
//@substitute: naming

static PyType_Slot __pyx_IterableCoroutineType_slots[] = {
    {Py_tp_dealloc, (void *)__Pyx_Coroutine_dealloc},
    {Py_am_await, (void *)&__Pyx_Coroutine_await},
    {Py_tp_traverse, (void *)__Pyx_Coroutine_traverse},
    {Py_tp_iter, (void *)__Pyx_Coroutine_await},
    {Py_tp_iternext, (void *)__Pyx_Generator_Next},
    {Py_tp_methods, (void *)__pyx_Coroutine_methods},
    {Py_tp_members, (void *)__pyx_Coroutine_memberlist},
    {Py_tp_getset, (void *)__pyx_Coroutine_getsets},
    {Py_tp_getattro, (void *) PyObject_GenericGetAttr},
#if CYTHON_USE_TP_FINALIZE
    {Py_tp_finalize, (void *)__Pyx_Coroutine_del},
#endif
#if __PYX_HAS_PY_AM_SEND == 1
    {Py_am_send, (void *)__Pyx_Coroutine_AmSend},
#endif
    {0, 0},
};

static PyType_Spec __pyx_IterableCoroutineType_spec = {
    __PYX_TYPE_MODULE_PREFIX "iterable_coroutine",
    sizeof(__pyx_CoroutineObject),
    0,
#if PY_VERSION_HEX >= 0x030C0000 && !CYTHON_COMPILING_IN_LIMITED_API
    Py_TPFLAGS_MANAGED_WEAKREF |
#endif
    Py_TPFLAGS_IMMUTABLETYPE | Py_TPFLAGS_DISALLOW_INSTANTIATION |
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | __Pyx_TPFLAGS_HAVE_AM_SEND, /*tp_flags*/
    __pyx_IterableCoroutineType_slots
};


static int __pyx_IterableCoroutine_init(PyObject *module) {
    $modulestatetype_cname *mstate = __Pyx_PyModule_GetState(module);
    mstate->__pyx_IterableCoroutineType = __Pyx_FetchCommonTypeFromSpec(
        mstate->__pyx_CommonTypesMetaclassType, module, &__pyx_IterableCoroutineType_spec, NULL);
    if (unlikely(!mstate->__pyx_IterableCoroutineType))
        return -1;
#if __PYX_HAS_PY_AM_SEND == 2
    __Pyx_SetBackportTypeAmSend(mstate->__pyx_IterableCoroutineType, &__pyx_Coroutine_as_async, &__Pyx_Coroutine_AmSend);
#endif
    return 0;
}

//////////////////// Generator.module_state_decls ////////////////////

PyTypeObject *__pyx_GeneratorType;

//////////////////// Generator.module_state_traverse ///////////////////

Py_VISIT(traverse_module_state->__pyx_GeneratorType);

//////////////////// Generator.module_state_clear ///////////////////

Py_CLEAR(clear_module_state->__pyx_GeneratorType);

//////////////////// Generator.init //////////////////
//@substitute: naming

if (likely(__pyx_Generator_init($module_cname) == 0)); else

//////////////////// Generator ////////////////////
//@requires: CoroutineBase
//@substitute: naming

static PyMethodDef __pyx_Generator_methods[] = {
    {"send", (PyCFunction) __Pyx_Coroutine_Send, METH_O,
     PyDoc_STR("send(arg) -> send 'arg' into generator,\nreturn next yielded value or raise StopIteration.")},
    {"throw", (PyCFunction) __Pyx_Coroutine_Throw, METH_VARARGS,
     PyDoc_STR("throw(typ[,val[,tb]]) -> raise exception in generator,\nreturn next yielded value or raise StopIteration.")},
    {"close", (PyCFunction) __Pyx_Coroutine_Close_Method, METH_NOARGS,
     PyDoc_STR("close() -> raise GeneratorExit inside generator.")},
    {"__reduce_ex__", (PyCFunction) __Pyx_Coroutine_fail_reduce_ex, METH_O, 0},
    {"__reduce__", (PyCFunction) __Pyx_Coroutine_fail_reduce_ex, METH_NOARGS, 0},
    {0, 0, 0, 0}
};

static PyMemberDef __pyx_Generator_memberlist[] = {
    {"gi_yieldfrom", T_OBJECT, offsetof(__pyx_CoroutineObject, yieldfrom), READONLY,
     PyDoc_STR("object being iterated by 'yield from', or None")},
    {"gi_code", T_OBJECT, offsetof(__pyx_CoroutineObject, gi_code), READONLY, NULL},
    {"__module__", T_OBJECT, offsetof(__pyx_CoroutineObject, gi_modulename), 0, 0},
#if PY_VERSION_HEX < 0x030C0000 || CYTHON_COMPILING_IN_LIMITED_API
    {"__weaklistoffset__", T_PYSSIZET, offsetof(__pyx_CoroutineObject, gi_weakreflist), READONLY, 0},
#endif
    {0, 0, 0, 0, 0}
};

static PyGetSetDef __pyx_Generator_getsets[] = {
    {"__name__", (getter)__Pyx_Coroutine_get_name, (setter)__Pyx_Coroutine_set_name,
     PyDoc_STR("name of the generator"), 0},
    {"__qualname__", (getter)__Pyx_Coroutine_get_qualname, (setter)__Pyx_Coroutine_set_qualname,
     PyDoc_STR("qualified name of the generator"), 0},
    {"gi_frame", (getter)__Pyx_Coroutine_get_frame, NULL,
     PyDoc_STR("Frame of the generator"), 0},
    {"gi_running", __Pyx_Coroutine_get_is_running_getter, NULL, NULL, NULL},
    {0, 0, 0, 0, 0}
};

static PyType_Slot __pyx_GeneratorType_slots[] = {
    {Py_tp_dealloc, (void *)__Pyx_Coroutine_dealloc},
    {Py_tp_traverse, (void *)__Pyx_Coroutine_traverse},
    {Py_tp_iter, (void *)PyObject_SelfIter},
    {Py_tp_iternext, (void *)__Pyx_Generator_Next},
    {Py_tp_methods, (void *)__pyx_Generator_methods},
    {Py_tp_members, (void *)__pyx_Generator_memberlist},
    {Py_tp_getset, (void *)__pyx_Generator_getsets},
    {Py_tp_getattro, (void *) PyObject_GenericGetAttr},
#if CYTHON_USE_TP_FINALIZE
    {Py_tp_finalize, (void *)__Pyx_Coroutine_del},
#endif
#if __PYX_HAS_PY_AM_SEND == 1
    {Py_am_send, (void *)__Pyx_Coroutine_AmSend},
#endif
    {0, 0},
};

static PyType_Spec __pyx_GeneratorType_spec = {
    __PYX_TYPE_MODULE_PREFIX "generator",
    sizeof(__pyx_CoroutineObject),
    0,
#if PY_VERSION_HEX >= 0x030C0000 && !CYTHON_COMPILING_IN_LIMITED_API
    Py_TPFLAGS_MANAGED_WEAKREF |
#endif
    Py_TPFLAGS_IMMUTABLETYPE | Py_TPFLAGS_DISALLOW_INSTANTIATION |
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | __Pyx_TPFLAGS_HAVE_AM_SEND, /*tp_flags*/
    __pyx_GeneratorType_slots
};

#if __PYX_HAS_PY_AM_SEND == 2
static __Pyx_PyAsyncMethodsStruct __pyx_Generator_as_async;
#endif

static int __pyx_Generator_init(PyObject *module) {
    $modulestatetype_cname *mstate = __Pyx_PyModule_GetState(module);
    mstate->__pyx_GeneratorType = __Pyx_FetchCommonTypeFromSpec(
        mstate->__pyx_CommonTypesMetaclassType, module, &__pyx_GeneratorType_spec, NULL);
    if (unlikely(!mstate->__pyx_GeneratorType)) {
        return -1;
    }
#if __PYX_HAS_PY_AM_SEND == 2
    __Pyx_SetBackportTypeAmSend(mstate->__pyx_GeneratorType, &__pyx_Generator_as_async, &__Pyx_Coroutine_AmSend);
#endif
    return 0;
}

static PyObject *__Pyx_Generator_GetInlinedResult(PyObject *self) {
    __pyx_CoroutineObject *gen = (__pyx_CoroutineObject*) self;
    PyObject *retval = NULL;
    if (unlikely(__Pyx_Coroutine_test_and_set_is_running(gen))) {
        return __Pyx_Coroutine_AlreadyRunningError(gen);
    }
    __Pyx_PySendResult result = __Pyx_Coroutine_SendEx(gen, Py_None, &retval, 0);
    __Pyx_Coroutine_unset_is_running(gen);
    (void) result;
    assert (result == PYGEN_RETURN || result == PYGEN_ERROR);
    assert ((result == PYGEN_RETURN && retval != NULL) || (result == PYGEN_ERROR && retval == NULL));
    return retval;
}


/////////////// ReturnWithStopIteration.proto ///////////////

static CYTHON_INLINE void __Pyx_ReturnWithStopIteration(PyObject* value, int async, int iternext); /*proto*/

/////////////// ReturnWithStopIteration ///////////////
//@requires: Exceptions.c::PyErrFetchRestore
//@requires: Exceptions.c::PyThreadStateGet
//@requires: ObjectHandling.c::PyObjectCallOneArg
//@substitute: naming

// 1) Instantiating an exception just to pass back a value is costly.
// 2) CPython 3.12 cannot separate exception type and value
// 3) Passing a tuple as value into PyErr_SetObject() passes its items on as arguments.
// 4) Passing an exception as value will interpret it as an exception on unpacking and raise it (or unpack its value).
// 5) If there is currently an exception being handled, we need to chain it.
// 6) CPython < 3.14a1 does not implement vectorcalls for exceptions.

static void __Pyx__ReturnWithStopIteration(PyObject* value, int async); /*proto*/

static CYTHON_INLINE void __Pyx_ReturnWithStopIteration(PyObject* value, int async, int iternext) {
    if (value == Py_None) {
        if (async || !iternext)
            PyErr_SetNone(async ? PyExc_StopAsyncIteration : PyExc_StopIteration);
        // for regular iternext, then we don't have to return StopIteration and can just return NULL
        return;
    }
    __Pyx__ReturnWithStopIteration(value, async);
}

static void __Pyx__ReturnWithStopIteration(PyObject* value, int async) {
#if CYTHON_COMPILING_IN_CPYTHON
    __Pyx_PyThreadState_declare
#endif
    PyObject *exc;
    PyObject *exc_type = async ? PyExc_StopAsyncIteration : PyExc_StopIteration;
#if CYTHON_COMPILING_IN_CPYTHON
    if ((PY_VERSION_HEX >= (0x030C00A6)) || unlikely(PyTuple_Check(value) || PyExceptionInstance_Check(value))) {
        if (PY_VERSION_HEX >= (0x030e00A1)) {
            exc = __Pyx_PyObject_CallOneArg(exc_type, value);
        } else {
            // Before Py3.14a1, CPython doesn't implement vectorcall for exceptions.
            PyObject *args_tuple = PyTuple_New(1);
            if (unlikely(!args_tuple)) return;
            Py_INCREF(value);
            PyTuple_SET_ITEM(args_tuple, 0, value);
            exc = PyObject_Call(exc_type, args_tuple, NULL);
            Py_DECREF(args_tuple);
        }
        if (unlikely(!exc)) return;
    } else {
        // it's safe to avoid instantiating the exception
        Py_INCREF(value);
        exc = value;
    }
    #if CYTHON_FAST_THREAD_STATE
    __Pyx_PyThreadState_assign
    #if CYTHON_USE_EXC_INFO_STACK
    if (!$local_tstate_cname->exc_info->exc_value)
    #else
    if (!$local_tstate_cname->exc_type)
    #endif
    {
        // no chaining needed => avoid the overhead in PyErr_SetObject()
        Py_INCREF(exc_type);
        __Pyx_ErrRestore(exc_type, exc, NULL);
        return;
    }
    #endif
#else
    exc = __Pyx_PyObject_CallOneArg(exc_type, value);
    if (unlikely(!exc)) return;
#endif
    PyErr_SetObject(exc_type, exc);
    Py_DECREF(exc);
}

//////////////////// Coro_CheckExact.proto ////////////////

#if CYTHON_COMPILING_IN_LIMITED_API
static int __Pyx_PyCoro_CheckExact(PyObject *o); /* proto */
#else
#define __Pyx_PyCoro_CheckExact PyCoro_CheckExact
#endif

/////////////////// Coro_CheckExact.module_state_decls //////////

#if CYTHON_COMPILING_IN_LIMITED_API
PyObject *__Pyx_CachedCoroType;
#endif

////////////////// Coro_CheckExact.init ////////////////

#if CYTHON_COMPILING_IN_LIMITED_API
{
    PyObject *typesModule=NULL;
    typesModule = PyImport_ImportModule("types");
    if (typesModule) {
        CGLOBAL(__Pyx_CachedCoroType) = PyObject_GetAttrString(typesModule, "CoroutineType");
        Py_DECREF(typesModule);
    }
} // error handling follows
#endif

/////////////// Coro_CheckExact.cleanup ////////////////

#if CYTHON_COMPILING_IN_LIMITED_API
Py_CLEAR(CGLOBAL(__Pyx_CachedCoroType));
#endif

/////////////////// Coro_CheckExact /////////////////////

#if CYTHON_COMPILING_IN_LIMITED_API
static int __Pyx_PyCoro_CheckExact(PyObject *o) {
    return (PyObject*)Py_TYPE(o) == CGLOBAL(__Pyx_CachedCoroType);
}
#endif
