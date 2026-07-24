// Exception raising code
//
// Exceptions are raised by __Pyx_Raise() and stored as plain
// type/value/tb in PyThreadState->curexc_*.  When being caught by an
// 'except' statement, curexc_* is moved over to exc_* by
// __Pyx_GetException()


/////////////// AssertionsEnabled.init ///////////////

if (likely(__Pyx_init_assertions_enabled() == 0)); else
// error propagation code is appended automatically

/////////////// AssertionsEnabled.proto ///////////////

#if CYTHON_COMPILING_IN_LIMITED_API  ||  PY_VERSION_HEX >= 0x030C0000
  // Py_OptimizeFlag is deprecated in Py3.12+ and not available in the Limited API.
  static int __pyx_assertions_enabled_flag;
  #define __pyx_assertions_enabled() (__pyx_assertions_enabled_flag)

  #if __clang__ || __GNUC__
  // "Assertions enabled" may be written multiple times when using subinterpreters.
  // However, it should always be written to the same value to isn't a "real" race.
  __attribute__((no_sanitize("thread")))
  #endif
  static int __Pyx_init_assertions_enabled(void) {
    PyObject *builtins, *debug, *debug_str;
    int flag;
    builtins = PyEval_GetBuiltins();
    if (!builtins) goto bad;
    // Not using PYIDENT() here because we probably don't need the string more than this once.
    debug_str = PyUnicode_FromStringAndSize("__debug__", 9);
    if (!debug_str) goto bad;
    debug = PyObject_GetItem(builtins, debug_str);
    Py_DECREF(debug_str);
    if (!debug) goto bad;
    flag = PyObject_IsTrue(debug);
    Py_DECREF(debug);
    if (flag == -1) goto bad;
    __pyx_assertions_enabled_flag = flag;
    return 0;
  bad:
    __pyx_assertions_enabled_flag = 1;
    // We (rarely) may not have an exception set, but the calling code will call PyErr_Occurred() either way.
    return -1;
  }
#else
  #define __Pyx_init_assertions_enabled()  (0)
  #define __pyx_assertions_enabled()  (!Py_OptimizeFlag)
#endif


/////////////// ErrOccurredWithGIL.proto ///////////////
static CYTHON_INLINE int __Pyx_ErrOccurredWithGIL(void); /* proto */

/////////////// ErrOccurredWithGIL ///////////////
static CYTHON_INLINE int __Pyx_ErrOccurredWithGIL(void) {
  int err;
  PyGILState_STATE _save = PyGILState_Ensure();
  err = !!PyErr_Occurred();
  PyGILState_Release(_save);
  return err;
}


/////////////// PyThreadStateGet.proto ///////////////
//@substitute: naming

#if CYTHON_FAST_THREAD_STATE
#define __Pyx_PyThreadState_declare  PyThreadState *$local_tstate_cname;
#define __Pyx_PyThreadState_assign  $local_tstate_cname = __Pyx_PyThreadState_Current;
#if PY_VERSION_HEX >= 0x030C00A6
#define __Pyx_PyErr_Occurred()  ($local_tstate_cname->current_exception != NULL)
#define __Pyx_PyErr_CurrentExceptionType()  ($local_tstate_cname->current_exception ? (PyObject*) Py_TYPE($local_tstate_cname->current_exception) : (PyObject*) NULL)
#else
#define __Pyx_PyErr_Occurred()  ($local_tstate_cname->curexc_type != NULL)
#define __Pyx_PyErr_CurrentExceptionType()  ($local_tstate_cname->curexc_type)
#endif
#else
// !CYTHON_FAST_THREAD_STATE
#define __Pyx_PyThreadState_declare
#define __Pyx_PyThreadState_assign
#define __Pyx_PyErr_Occurred()  (PyErr_Occurred() != NULL)
#define __Pyx_PyErr_CurrentExceptionType()  PyErr_Occurred()
#endif


/////////////// PyErrExceptionMatches.proto ///////////////
//@substitute: naming

#if CYTHON_FAST_THREAD_STATE
#define __Pyx_PyErr_ExceptionMatches(err) __Pyx_PyErr_ExceptionMatchesInState($local_tstate_cname, err)
static CYTHON_INLINE int __Pyx_PyErr_ExceptionMatchesInState(PyThreadState* tstate, PyObject* err);
#else
#define __Pyx_PyErr_ExceptionMatches(err)  PyErr_ExceptionMatches(err)
#endif

/////////////// PyErrExceptionMatches ///////////////

#if CYTHON_FAST_THREAD_STATE
static int __Pyx_PyErr_ExceptionMatchesTuple(PyObject *exc_type, PyObject *tuple) {
    Py_ssize_t i, n;
    n = PyTuple_GET_SIZE(tuple);
    // the tighter subtype checking in Py3 allows faster out-of-order comparison
    for (i=0; i<n; i++) {
        if (exc_type == PyTuple_GET_ITEM(tuple, i)) return 1;
    }
    for (i=0; i<n; i++) {
        if (__Pyx_PyErr_GivenExceptionMatches(exc_type, PyTuple_GET_ITEM(tuple, i))) return 1;
    }
    return 0;
}

static CYTHON_INLINE int __Pyx_PyErr_ExceptionMatchesInState(PyThreadState* tstate, PyObject* err) {
    int result;
    PyObject *exc_type;
#if PY_VERSION_HEX >= 0x030C00A6
    PyObject *current_exception = tstate->current_exception;
    if (unlikely(!current_exception)) return 0;
    exc_type = (PyObject*) Py_TYPE(current_exception);
    if (exc_type == err) return 1;
#else
    exc_type = tstate->curexc_type;
    if (exc_type == err) return 1;
    if (unlikely(!exc_type)) return 0;
#endif
    #if CYTHON_AVOID_BORROWED_REFS
    Py_INCREF(exc_type);
    #endif
    if (unlikely(PyTuple_Check(err))) {
        result = __Pyx_PyErr_ExceptionMatchesTuple(exc_type, err);
    } else {
        result = __Pyx_PyErr_GivenExceptionMatches(exc_type, err);
    }
    #if CYTHON_AVOID_BORROWED_REFS
    Py_DECREF(exc_type);
    #endif
    return result;
}
#endif

/////////////// PyErrFetchRestore.proto ///////////////
//@substitute: naming
//@requires: PyThreadStateGet

#if CYTHON_FAST_THREAD_STATE
#define __Pyx_PyErr_Clear() __Pyx_ErrRestore(NULL, NULL, NULL)
#define __Pyx_ErrRestoreWithState(type, value, tb)  __Pyx_ErrRestoreInState(PyThreadState_GET(), type, value, tb)
#define __Pyx_ErrFetchWithState(type, value, tb)    __Pyx_ErrFetchInState(PyThreadState_GET(), type, value, tb)
#define __Pyx_ErrRestore(type, value, tb)  __Pyx_ErrRestoreInState($local_tstate_cname, type, value, tb)
#define __Pyx_ErrFetch(type, value, tb)    __Pyx_ErrFetchInState($local_tstate_cname, type, value, tb)
static CYTHON_INLINE void __Pyx_ErrRestoreInState(PyThreadState *tstate, PyObject *type, PyObject *value, PyObject *tb); /*proto*/
static CYTHON_INLINE void __Pyx_ErrFetchInState(PyThreadState *tstate, PyObject **type, PyObject **value, PyObject **tb); /*proto*/

#if CYTHON_COMPILING_IN_CPYTHON && PY_VERSION_HEX < 0x030C00A6
#define __Pyx_PyErr_SetNone(exc) (Py_INCREF(exc), __Pyx_ErrRestore((exc), NULL, NULL))
#else
#define __Pyx_PyErr_SetNone(exc) PyErr_SetNone(exc)
#endif

#else
#define __Pyx_PyErr_Clear() PyErr_Clear()
#define __Pyx_PyErr_SetNone(exc) PyErr_SetNone(exc)
#define __Pyx_ErrRestoreWithState(type, value, tb)  PyErr_Restore(type, value, tb)
#define __Pyx_ErrFetchWithState(type, value, tb)  PyErr_Fetch(type, value, tb)
#define __Pyx_ErrRestoreInState(tstate, type, value, tb)  PyErr_Restore(type, value, tb)
#define __Pyx_ErrFetchInState(tstate, type, value, tb)  PyErr_Fetch(type, value, tb)
#define __Pyx_ErrRestore(type, value, tb)  PyErr_Restore(type, value, tb)
#define __Pyx_ErrFetch(type, value, tb)  PyErr_Fetch(type, value, tb)
#endif

/////////////// PyErrFetchRestore ///////////////

#if CYTHON_FAST_THREAD_STATE
static CYTHON_INLINE void __Pyx_ErrRestoreInState(PyThreadState *tstate, PyObject *type, PyObject *value, PyObject *tb) {
#if PY_VERSION_HEX >= 0x030C00A6
    PyObject *tmp_value;
    assert(type == NULL || (value != NULL && type == (PyObject*) Py_TYPE(value)));
    if (value) {
        #if CYTHON_COMPILING_IN_CPYTHON
        if (unlikely(((PyBaseExceptionObject*) value)->traceback != tb))
        #endif
            // If this fails, we may lose the traceback but still set the expected exception below.
            PyException_SetTraceback(value, tb);
    }
    tmp_value = tstate->current_exception;
    tstate->current_exception = value;
    Py_XDECREF(tmp_value);
    Py_XDECREF(type);
    Py_XDECREF(tb);
#else
    PyObject *tmp_type, *tmp_value, *tmp_tb;
    tmp_type = tstate->curexc_type;
    tmp_value = tstate->curexc_value;
    tmp_tb = tstate->curexc_traceback;
    tstate->curexc_type = type;
    tstate->curexc_value = value;
    tstate->curexc_traceback = tb;
    Py_XDECREF(tmp_type);
    Py_XDECREF(tmp_value);
    Py_XDECREF(tmp_tb);
#endif
}

static CYTHON_INLINE void __Pyx_ErrFetchInState(PyThreadState *tstate, PyObject **type, PyObject **value, PyObject **tb) {
#if PY_VERSION_HEX >= 0x030C00A6
    PyObject* exc_value;
    exc_value = tstate->current_exception;
    tstate->current_exception = 0;
    *value = exc_value;
    *type = NULL;
    *tb = NULL;
    if (exc_value) {
        *type = (PyObject*) Py_TYPE(exc_value);
        Py_INCREF(*type);
        #if CYTHON_COMPILING_IN_CPYTHON
        *tb = ((PyBaseExceptionObject*) exc_value)->traceback;
        Py_XINCREF(*tb);
        #else
        *tb = PyException_GetTraceback(exc_value);
        #endif
    }
#else
    *type = tstate->curexc_type;
    *value = tstate->curexc_value;
    *tb = tstate->curexc_traceback;
    tstate->curexc_type = 0;
    tstate->curexc_value = 0;
    tstate->curexc_traceback = 0;
#endif
}
#endif

/////////////// RaiseException.export ///////////////

static void __Pyx_Raise(PyObject *type, PyObject *value, PyObject *tb, PyObject *cause); /*proto*/

/////////////// RaiseException ///////////////
//@requires: PyErrFetchRestore
//@requires: PyThreadStateGet

// The following function is based on do_raise() from ceval.c.

static void __Pyx_Raise(PyObject *type, PyObject *value, PyObject *tb, PyObject *cause) {
    PyObject* owned_instance = NULL;
    if (tb == Py_None) {
        tb = 0;
    } else if (tb && !PyTraceBack_Check(tb)) {
        PyErr_SetString(PyExc_TypeError,
            "raise: arg 3 must be a traceback or None");
        goto bad;
    }
    if (value == Py_None)
        value = 0;

    if (PyExceptionInstance_Check(type)) {
        if (value) {
            PyErr_SetString(PyExc_TypeError,
                "instance exception may not have a separate value");
            goto bad;
        }
        value = type;
        type = (PyObject*) Py_TYPE(value);
    } else if (PyExceptionClass_Check(type)) {
        // make sure value is an exception instance of type
        PyObject *instance_class = NULL;
        if (value && PyExceptionInstance_Check(value)) {
            instance_class = (PyObject*) Py_TYPE(value);
            if (instance_class != type) {
                int is_subclass = PyObject_IsSubclass(instance_class, type);
                if (!is_subclass) {
                    instance_class = NULL;
                } else if (unlikely(is_subclass == -1)) {
                    // error on subclass test
                    goto bad;
                } else {
                    // believe the instance
                    type = instance_class;
                }
            }
        }
        if (!instance_class) {
            // instantiate the type now (we don't know when and how it will be caught)
            // assuming that 'value' is an argument to the type's constructor
            // not using PyErr_NormalizeException() to avoid ref-counting problems
            PyObject *args;
            if (!value)
                args = PyTuple_New(0);
            else if (PyTuple_Check(value)) {
                Py_INCREF(value);
                args = value;
            } else
                args = PyTuple_Pack(1, value);
            if (!args)
                goto bad;
            owned_instance = PyObject_Call(type, args, NULL);
            Py_DECREF(args);
            if (!owned_instance)
                goto bad;
            value = owned_instance;
            if (!PyExceptionInstance_Check(value)) {
                PyErr_Format(PyExc_TypeError,
                             "calling %R should have returned an instance of "
                             "BaseException, not %R",
                             type, Py_TYPE(value));
                goto bad;
            }
        }
    } else {
        PyErr_SetString(PyExc_TypeError,
            "raise: exception class must be a subclass of BaseException");
        goto bad;
    }

    if (cause) {
        PyObject *fixed_cause;
        if (cause == Py_None) {
            // raise ... from None
            fixed_cause = NULL;
        } else if (PyExceptionClass_Check(cause)) {
            fixed_cause = PyObject_CallObject(cause, NULL);
            if (fixed_cause == NULL)
                goto bad;
        } else if (PyExceptionInstance_Check(cause)) {
            fixed_cause = cause;
            Py_INCREF(fixed_cause);
        } else {
            PyErr_SetString(PyExc_TypeError,
                            "exception causes must derive from "
                            "BaseException");
            goto bad;
        }
        PyException_SetCause(value, fixed_cause);
    }

    PyErr_SetObject(type, value);

    if (tb) {
#if PY_VERSION_HEX >= 0x030C00A6
        // If this fails, we just get a different exception, so ignore the return value.
        PyException_SetTraceback(value, tb);
#elif CYTHON_FAST_THREAD_STATE
        PyThreadState *tstate = __Pyx_PyThreadState_Current;
        PyObject* tmp_tb = tstate->curexc_traceback;
        if (tb != tmp_tb) {
            Py_INCREF(tb);
            tstate->curexc_traceback = tb;
            Py_XDECREF(tmp_tb);
        }
#else
        PyObject *tmp_type, *tmp_value, *tmp_tb;
        PyErr_Fetch(&tmp_type, &tmp_value, &tmp_tb);
        Py_INCREF(tb);
        PyErr_Restore(tmp_type, tmp_value, tb);
        Py_XDECREF(tmp_tb);
#endif
    }

bad:
    Py_XDECREF(owned_instance);
    return;
}


/////////////// GetTopmostException.proto ///////////////

#if CYTHON_USE_EXC_INFO_STACK && CYTHON_FAST_THREAD_STATE
static _PyErr_StackItem * __Pyx_PyErr_GetTopmostException(PyThreadState *tstate);
#endif

/////////////// GetTopmostException ///////////////

#if CYTHON_USE_EXC_INFO_STACK && CYTHON_FAST_THREAD_STATE
// Copied from errors.c in CPython.
static _PyErr_StackItem *
__Pyx_PyErr_GetTopmostException(PyThreadState *tstate)
{
    _PyErr_StackItem *exc_info = tstate->exc_info;
    while ((exc_info->exc_value == NULL || exc_info->exc_value == Py_None) &&
           exc_info->previous_item != NULL)
    {
        exc_info = exc_info->previous_item;
    }
    return exc_info;
}
#endif


/////////////// GetException.proto ///////////////
//@substitute: naming
//@requires: PyThreadStateGet

#if CYTHON_FAST_THREAD_STATE
#define __Pyx_GetException(type, value, tb)  __Pyx__GetException($local_tstate_cname, type, value, tb)
static int __Pyx__GetException(PyThreadState *tstate, PyObject **type, PyObject **value, PyObject **tb); /*proto*/
#else
static int __Pyx_GetException(PyObject **type, PyObject **value, PyObject **tb); /*proto*/
#endif

/////////////// GetException ///////////////

#if CYTHON_FAST_THREAD_STATE
static int __Pyx__GetException(PyThreadState *tstate, PyObject **type, PyObject **value, PyObject **tb)
#else
static int __Pyx_GetException(PyObject **type, PyObject **value, PyObject **tb)
#endif
{
    PyObject *local_type = NULL, *local_value, *local_tb = NULL;
#if CYTHON_FAST_THREAD_STATE
    PyObject *tmp_type, *tmp_value, *tmp_tb;
  #if PY_VERSION_HEX >= 0x030C0000
    local_value = tstate->current_exception;
    tstate->current_exception = 0;
  #else
    local_type = tstate->curexc_type;
    local_value = tstate->curexc_value;
    local_tb = tstate->curexc_traceback;
    tstate->curexc_type = 0;
    tstate->curexc_value = 0;
    tstate->curexc_traceback = 0;
  #endif
#elif __PYX_LIMITED_VERSION_HEX > 0x030C0000
    local_value = PyErr_GetRaisedException();
#else
    PyErr_Fetch(&local_type, &local_value, &local_tb);
#endif

#if __PYX_LIMITED_VERSION_HEX > 0x030C0000
    if (likely(local_value)) {
        local_type = (PyObject*) Py_TYPE(local_value);
        Py_INCREF(local_type);
        local_tb = PyException_GetTraceback(local_value);
    }
#else
    // Note that In Python 3.12+ exceptions are already normalized
    PyErr_NormalizeException(&local_type, &local_value, &local_tb);
#if CYTHON_FAST_THREAD_STATE
    if (unlikely(tstate->curexc_type))
#else
    if (unlikely(PyErr_Occurred()))
#endif
        goto bad;

    // Note that in Python 3.12 the traceback came directly from local_value anyway.
    if (local_tb) {
        if (unlikely(PyException_SetTraceback(local_value, local_tb) < 0))
            goto bad;
    }
#endif // __PYX_LIMITED_VERSION_HEX > 0x030C0000

    // traceback may be NULL for freshly raised exceptions
    Py_XINCREF(local_tb);
    // exception state may be empty in parallel loops (code-gen error where we don't generate a 
    // top-level PyGILState_Ensure surrounding the whole loop, and so releasing the GIL temporarily
    // wipes the whole thread state).
    Py_XINCREF(local_type);
    Py_XINCREF(local_value);
    *type = local_type;
    *value = local_value;
    *tb = local_tb;

#if CYTHON_FAST_THREAD_STATE
    #if CYTHON_USE_EXC_INFO_STACK
    {
        _PyErr_StackItem *exc_info = tstate->exc_info;
      #if PY_VERSION_HEX >= 0x030B00a4
        tmp_value = exc_info->exc_value;
        exc_info->exc_value = local_value;
        tmp_type = NULL;
        tmp_tb = NULL;
        Py_XDECREF(local_type);
        Py_XDECREF(local_tb);
      #else
        tmp_type = exc_info->exc_type;
        tmp_value = exc_info->exc_value;
        tmp_tb = exc_info->exc_traceback;
        exc_info->exc_type = local_type;
        exc_info->exc_value = local_value;
        exc_info->exc_traceback = local_tb;
      #endif
    }
    #else
    tmp_type = tstate->exc_type;
    tmp_value = tstate->exc_value;
    tmp_tb = tstate->exc_traceback;
    tstate->exc_type = local_type;
    tstate->exc_value = local_value;
    tstate->exc_traceback = local_tb;
    #endif
    // Make sure tstate is in a consistent state when we XDECREF
    // these objects (DECREF may run arbitrary code).
    Py_XDECREF(tmp_type);
    Py_XDECREF(tmp_value);
    Py_XDECREF(tmp_tb);
#elif __PYX_LIMITED_VERSION_HEX >= 0x030b0000
    PyErr_SetHandledException(local_value);
    Py_XDECREF(local_value);
    Py_XDECREF(local_type);
    Py_XDECREF(local_tb);
#else
    PyErr_SetExcInfo(local_type, local_value, local_tb);
#endif

    return 0;

#if __PYX_LIMITED_VERSION_HEX <= 0x030C0000
bad:
    *type = 0;
    *value = 0;
    *tb = 0;
    Py_XDECREF(local_type);
    Py_XDECREF(local_value);
    Py_XDECREF(local_tb);
    return -1;
#endif
}

/////////////// ReRaiseException.proto ///////////////

static CYTHON_INLINE void __Pyx_ReraiseException(void); /*proto*/

/////////////// ReRaiseException ///////////////
//@requires: GetTopmostException

static CYTHON_INLINE void __Pyx_ReraiseException(void) {
    PyObject *type = NULL, *value = NULL, *tb = NULL;
#if CYTHON_FAST_THREAD_STATE
    PyThreadState *tstate = PyThreadState_GET();
  #if CYTHON_USE_EXC_INFO_STACK
    _PyErr_StackItem *exc_info = __Pyx_PyErr_GetTopmostException(tstate);
    value = exc_info->exc_value;
    #if PY_VERSION_HEX >= 0x030B00a4
    if (unlikely(value == Py_None)) {
        value = NULL;
    } else if (value) {
        Py_INCREF(value);
        type = (PyObject*) Py_TYPE(value);
        Py_INCREF(type);
        tb = PyException_GetTraceback(value);
    }
    #else
    type = exc_info->exc_type;
    tb = exc_info->exc_traceback;
    Py_XINCREF(type);
    Py_XINCREF(value);
    Py_XINCREF(tb);
    #endif
  #else
    type = tstate->exc_type;
    value = tstate->exc_value;
    tb = tstate->exc_traceback;
    Py_XINCREF(type);
    Py_XINCREF(value);
    Py_XINCREF(tb);
  #endif
#else
    PyErr_GetExcInfo(&type, &value, &tb);
#endif
    if (unlikely(!type || type == Py_None)) {
        Py_XDECREF(type);
        Py_XDECREF(value);
        Py_XDECREF(tb);
        // message copied from Py3
        PyErr_SetString(PyExc_RuntimeError,
            "No active exception to reraise");
    } else {
        PyErr_Restore(type, value, tb);
    }
}

/////////////// SaveResetException.proto ///////////////
//@substitute: naming
//@requires: PyThreadStateGet

#if CYTHON_FAST_THREAD_STATE
#define __Pyx_ExceptionSave(type, value, tb)  __Pyx__ExceptionSave($local_tstate_cname, type, value, tb)
static CYTHON_INLINE void __Pyx__ExceptionSave(PyThreadState *tstate, PyObject **type, PyObject **value, PyObject **tb); /*proto*/
#define __Pyx_ExceptionReset(type, value, tb)  __Pyx__ExceptionReset($local_tstate_cname, type, value, tb)
static CYTHON_INLINE void __Pyx__ExceptionReset(PyThreadState *tstate, PyObject *type, PyObject *value, PyObject *tb); /*proto*/

#else

#define __Pyx_ExceptionSave(type, value, tb)   PyErr_GetExcInfo(type, value, tb)
#define __Pyx_ExceptionReset(type, value, tb)  PyErr_SetExcInfo(type, value, tb)
#endif

/////////////// SaveResetException ///////////////
//@requires: GetTopmostException

#if CYTHON_FAST_THREAD_STATE
static CYTHON_INLINE void __Pyx__ExceptionSave(PyThreadState *tstate, PyObject **type, PyObject **value, PyObject **tb) {
  #if CYTHON_USE_EXC_INFO_STACK && PY_VERSION_HEX >= 0x030B00a4
    _PyErr_StackItem *exc_info = __Pyx_PyErr_GetTopmostException(tstate);
    PyObject *exc_value = exc_info->exc_value;
    if (exc_value == NULL || exc_value == Py_None) {
        *value = NULL;
        *type = NULL;
        *tb = NULL;
    } else {
        *value = exc_value;
        Py_INCREF(*value);
        *type = (PyObject*) Py_TYPE(exc_value);
        Py_INCREF(*type);
        *tb = PyException_GetTraceback(exc_value);
    }
  #elif CYTHON_USE_EXC_INFO_STACK
    _PyErr_StackItem *exc_info = __Pyx_PyErr_GetTopmostException(tstate);
    *type = exc_info->exc_type;
    *value = exc_info->exc_value;
    *tb = exc_info->exc_traceback;
    Py_XINCREF(*type);
    Py_XINCREF(*value);
    Py_XINCREF(*tb);
  #else
    *type = tstate->exc_type;
    *value = tstate->exc_value;
    *tb = tstate->exc_traceback;
    Py_XINCREF(*type);
    Py_XINCREF(*value);
    Py_XINCREF(*tb);
  #endif
}

static CYTHON_INLINE void __Pyx__ExceptionReset(PyThreadState *tstate, PyObject *type, PyObject *value, PyObject *tb) {
  #if CYTHON_USE_EXC_INFO_STACK && PY_VERSION_HEX >= 0x030B00a4
    _PyErr_StackItem *exc_info = tstate->exc_info;
    PyObject *tmp_value = exc_info->exc_value;
    exc_info->exc_value = value;
    Py_XDECREF(tmp_value);
    // TODO: avoid passing these at all
    Py_XDECREF(type);
    Py_XDECREF(tb);
  #else
    PyObject *tmp_type, *tmp_value, *tmp_tb;
    #if CYTHON_USE_EXC_INFO_STACK
    _PyErr_StackItem *exc_info = tstate->exc_info;
    tmp_type = exc_info->exc_type;
    tmp_value = exc_info->exc_value;
    tmp_tb = exc_info->exc_traceback;
    exc_info->exc_type = type;
    exc_info->exc_value = value;
    exc_info->exc_traceback = tb;
    #else
    tmp_type = tstate->exc_type;
    tmp_value = tstate->exc_value;
    tmp_tb = tstate->exc_traceback;
    tstate->exc_type = type;
    tstate->exc_value = value;
    tstate->exc_traceback = tb;
    #endif
    Py_XDECREF(tmp_type);
    Py_XDECREF(tmp_value);
    Py_XDECREF(tmp_tb);
  #endif
}
#endif

/////////////// SwapException.proto ///////////////
//@substitute: naming
//@requires: PyThreadStateGet

#if CYTHON_FAST_THREAD_STATE
#define __Pyx_ExceptionSwap(type, value, tb)  __Pyx__ExceptionSwap($local_tstate_cname, type, value, tb)
static CYTHON_INLINE void __Pyx__ExceptionSwap(PyThreadState *tstate, PyObject **type, PyObject **value, PyObject **tb); /*proto*/
#else
static CYTHON_INLINE void __Pyx_ExceptionSwap(PyObject **type, PyObject **value, PyObject **tb); /*proto*/
#endif

/////////////// SwapException ///////////////

#if CYTHON_FAST_THREAD_STATE
static CYTHON_INLINE void __Pyx__ExceptionSwap(PyThreadState *tstate, PyObject **type, PyObject **value, PyObject **tb) {
    PyObject *tmp_type, *tmp_value, *tmp_tb;
  #if CYTHON_USE_EXC_INFO_STACK && PY_VERSION_HEX >= 0x030B00a4
    _PyErr_StackItem *exc_info = tstate->exc_info;
    tmp_value = exc_info->exc_value;
    exc_info->exc_value = *value;
    if (tmp_value == NULL || tmp_value == Py_None) {
        Py_XDECREF(tmp_value);
        tmp_value = NULL;
        tmp_type = NULL;
        tmp_tb = NULL;
    } else {
        // TODO: avoid swapping these at all
        tmp_type = (PyObject*) Py_TYPE(tmp_value);
        Py_INCREF(tmp_type);
        #if CYTHON_COMPILING_IN_CPYTHON
        tmp_tb = ((PyBaseExceptionObject*) tmp_value)->traceback;
        Py_XINCREF(tmp_tb);
        #else
        tmp_tb = PyException_GetTraceback(tmp_value);
        #endif
    }
  #elif CYTHON_USE_EXC_INFO_STACK
    _PyErr_StackItem *exc_info = tstate->exc_info;
    tmp_type = exc_info->exc_type;
    tmp_value = exc_info->exc_value;
    tmp_tb = exc_info->exc_traceback;

    exc_info->exc_type = *type;
    exc_info->exc_value = *value;
    exc_info->exc_traceback = *tb;
  #else
    tmp_type = tstate->exc_type;
    tmp_value = tstate->exc_value;
    tmp_tb = tstate->exc_traceback;

    tstate->exc_type = *type;
    tstate->exc_value = *value;
    tstate->exc_traceback = *tb;
  #endif

    *type = tmp_type;
    *value = tmp_value;
    *tb = tmp_tb;
}

#else

static CYTHON_INLINE void __Pyx_ExceptionSwap(PyObject **type, PyObject **value, PyObject **tb) {
    PyObject *tmp_type, *tmp_value, *tmp_tb;
    PyErr_GetExcInfo(&tmp_type, &tmp_value, &tmp_tb);
    PyErr_SetExcInfo(*type, *value, *tb);
    *type = tmp_type;
    *value = tmp_value;
    *tb = tmp_tb;
}
#endif

/////////////// WriteUnraisableException.proto ///////////////

static void __Pyx_WriteUnraisable(const char *name, int clineno,
                                  int lineno, const char *filename,
                                  int full_traceback, int nogil); /*proto*/

/////////////// WriteUnraisableException ///////////////
//@requires: PyErrFetchRestore
//@requires: PyThreadStateGet

static void __Pyx_WriteUnraisable(const char *name, int clineno,
                                  int lineno, const char *filename,
                                  int full_traceback, int nogil) {
    PyObject *old_exc, *old_val, *old_tb;
    PyObject *ctx;
    __Pyx_PyThreadState_declare
    PyGILState_STATE state;
    if (nogil)
        state = PyGILState_Ensure();
    /* arbitrary, to suppress warning */
    else state = (PyGILState_STATE)0;
    CYTHON_UNUSED_VAR(clineno);
    CYTHON_UNUSED_VAR(lineno);
    CYTHON_UNUSED_VAR(filename);
    CYTHON_MAYBE_UNUSED_VAR(nogil);

    __Pyx_PyThreadState_assign
    __Pyx_ErrFetch(&old_exc, &old_val, &old_tb);
    if (full_traceback) {
        Py_XINCREF(old_exc);
        Py_XINCREF(old_val);
        Py_XINCREF(old_tb);
        __Pyx_ErrRestore(old_exc, old_val, old_tb);
        PyErr_PrintEx(0);
    }
    ctx = PyUnicode_FromString(name);
    __Pyx_ErrRestore(old_exc, old_val, old_tb);
    if (!ctx) {
        PyErr_WriteUnraisable(Py_None);
    } else {
        PyErr_WriteUnraisable(ctx);
        Py_DECREF(ctx);
    }
    if (nogil)
        PyGILState_Release(state);
}

/////////////// CLineInTraceback.proto ///////////////

#if CYTHON_CLINE_IN_TRACEBACK && CYTHON_CLINE_IN_TRACEBACK_RUNTIME
static int __Pyx_CLineForTraceback(PyThreadState *tstate, int c_line);/*proto*/
#else
#define __Pyx_CLineForTraceback(tstate, c_line)  (((CYTHON_CLINE_IN_TRACEBACK)) ? c_line : 0)
#endif

/////////////// CLineInTraceback ///////////////
//@requires: ObjectHandling.c::PyObjectGetAttrStrNoError
//@requires: ObjectHandling.c::PyDictVersioning
//@requires: PyErrFetchRestore
//@requires: Builtins.c::dict_setdefault

#if CYTHON_CLINE_IN_TRACEBACK && CYTHON_CLINE_IN_TRACEBACK_RUNTIME
#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000
// On earlier version of the Limited API we're less flexible and assume it's
// definitely a module. Which will break some odd user monkey-patching...
#define __Pyx_PyProbablyModule_GetDict(o) __Pyx_XNewRef(PyModule_GetDict(o))
#elif !CYTHON_COMPILING_IN_CPYTHON || CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
// On freethreading, use PyObject_GenericGetDict to avoid having to think too hard
// about atomics and _PyObject_GetDictPtr.
#define __Pyx_PyProbablyModule_GetDict(o) PyObject_GenericGetDict(o, NULL);
#else
PyObject* __Pyx_PyProbablyModule_GetDict(PyObject *o) {
    PyObject **dict_ptr = _PyObject_GetDictPtr(o);
    return dict_ptr ? __Pyx_XNewRef(*dict_ptr) : NULL;
}
#endif


static int __Pyx_CLineForTraceback(PyThreadState *tstate, int c_line) {
    PyObject *use_cline = NULL;
    PyObject *ptype, *pvalue, *ptraceback;
    PyObject *cython_runtime_dict;

    CYTHON_MAYBE_UNUSED_VAR(tstate);

    if (unlikely(!NAMED_CGLOBAL(cython_runtime_cname))) {
        // Very early error where the runtime module is not set up yet.
        return c_line;
    }

    __Pyx_ErrFetchInState(tstate, &ptype, &pvalue, &ptraceback);

    cython_runtime_dict = __Pyx_PyProbablyModule_GetDict(NAMED_CGLOBAL(cython_runtime_cname));
    if (likely(cython_runtime_dict)) {
        __PYX_PY_DICT_LOOKUP_IF_MODIFIED(
            use_cline, cython_runtime_dict,
            __Pyx_PyDict_SetDefault(cython_runtime_dict, PYIDENT("cline_in_traceback"), Py_False))
    }
    if (use_cline == NULL || use_cline == Py_False || (use_cline != Py_True && PyObject_Not(use_cline) != 0)) {
        c_line = 0;
    }
    Py_XDECREF(use_cline);
    Py_XDECREF(cython_runtime_dict);
    __Pyx_ErrRestoreInState(tstate, ptype, pvalue, ptraceback);
    return c_line;
}
#endif

/////////////// AddTraceback.proto ///////////////

static void __Pyx_AddTraceback(const char *funcname, int c_line,
                               int py_line, const char *filename); /*proto*/

/////////////// AddTraceback ///////////////
//@requires: ModuleSetupCode.c::CodeObjectCache
//@requires: CLineInTraceback
//@substitute: naming

#include "compile.h"
#include "frameobject.h"
#include "traceback.h"
#if PY_VERSION_HEX >= 0x030b00a6 && !CYTHON_COMPILING_IN_LIMITED_API && !defined(PYPY_VERSION)
  #ifndef Py_BUILD_CORE
    #define Py_BUILD_CORE 1
  #endif
  #include "internal/pycore_frame.h"
#endif

#if CYTHON_COMPILING_IN_LIMITED_API
static PyObject *__Pyx_PyCode_Replace_For_AddTraceback(PyObject *code, PyObject *scratch_dict,
                                                       PyObject *firstlineno, PyObject *name) {
    PyObject *replace = NULL;
    if (unlikely(PyDict_SetItemString(scratch_dict, "co_firstlineno", firstlineno))) return NULL;
    if (unlikely(PyDict_SetItemString(scratch_dict, "co_name", name))) return NULL;

    replace = PyObject_GetAttrString(code, "replace");
    if (likely(replace)) {
        PyObject *result = PyObject_Call(replace, EMPTY(tuple), scratch_dict);
        Py_DECREF(replace);
        return result;
    }
    PyErr_Clear();

    return NULL;
}

static void __Pyx_AddTraceback(const char *funcname, int c_line,
                               int py_line, const char *filename) {
    PyObject *code_object = NULL, *py_py_line = NULL, *py_funcname = NULL, *dict = NULL;
    PyObject *replace = NULL, *getframe = NULL, *frame = NULL;
    PyObject *exc_type, *exc_value, *exc_traceback;
    int success = 0;
    if (c_line) {
        c_line = __Pyx_CLineForTraceback(__Pyx_PyThreadState_Current, c_line);
    }

    // DW - this is a horrendous hack, but I'm quite proud of it. Essentially
    // we need to generate a frame with the right line number/filename/funcname.
    // We do this by compiling a small bit of code that uses sys._getframe to get a
    // frame, and then customizing the details of the code to match.
    // We then run the code object and use the generated frame to set the traceback.

    PyErr_Fetch(&exc_type, &exc_value, &exc_traceback);

    code_object = $global_code_object_cache_find(c_line ? -c_line : py_line);
    if (!code_object) {
        code_object = Py_CompileString("_getframe()", filename, Py_eval_input);
        if (unlikely(!code_object)) goto bad;
        py_py_line = PyLong_FromLong(py_line);
        if (unlikely(!py_py_line)) goto bad;
        if (c_line) {
            py_funcname = PyUnicode_FromFormat( "%s (%s:%d)", funcname, $cfilenm_cname, c_line);
        } else {
            py_funcname = PyUnicode_FromString(funcname);
        }
        if (unlikely(!py_funcname)) goto bad;
        dict = PyDict_New();
        if (unlikely(!dict)) goto bad;
        {
            PyObject *old_code_object = code_object;
            code_object = __Pyx_PyCode_Replace_For_AddTraceback(code_object, dict, py_py_line, py_funcname);
            Py_DECREF(old_code_object);
        }
        if (unlikely(!code_object)) goto bad;

        $global_code_object_cache_insert(c_line ? -c_line : py_line, code_object);
    } else {
        // The frame part still expects a dict
        dict = PyDict_New();
    }

    // Note that getframe is borrowed
    getframe = PySys_GetObject("_getframe");
    if (unlikely(!getframe)) goto bad;
    // reuse dict as globals (nothing conflicts, and it saves an allocation)
    if (unlikely(PyDict_SetItemString(dict, "_getframe", getframe))) goto bad;

    frame = PyEval_EvalCode(code_object, dict, dict);
    if (unlikely(!frame) || frame == Py_None) goto bad;
    success = 1;

  bad:
    PyErr_Restore(exc_type, exc_value, exc_traceback);
    Py_XDECREF(code_object);
    Py_XDECREF(py_py_line);
    Py_XDECREF(py_funcname);
    Py_XDECREF(dict);
    Py_XDECREF(replace);

    if (success) {
        // Unfortunately an easy way to check the type of frame isn't in the
        // limited API. The check against None should cover the most
        // likely wrong answer though.
        PyTraceBack_Here(
            // Python < 0x03090000 didn't expose PyFrameObject
            // but they all expose struct _frame as an opaque type
            (struct _frame*)frame);
    }

    Py_XDECREF(frame);
}
#else
static PyCodeObject* __Pyx_CreateCodeObjectForTraceback(
            const char *funcname, int c_line,
            int py_line, const char *filename) {
    PyCodeObject *py_code = NULL;
    PyObject *py_funcname = NULL;

    if (c_line) {
        py_funcname = PyUnicode_FromFormat( "%s (%s:%d)", funcname, $cfilenm_cname, c_line);
        if (!py_funcname) goto bad;
        funcname = PyUnicode_AsUTF8(py_funcname);
        if (!funcname) goto bad;
    }
    py_code = PyCode_NewEmpty(filename, funcname, py_line);
    Py_XDECREF(py_funcname);  /* XDECREF since it's only set on Py3 if cline */
    return py_code;
bad:
    Py_XDECREF(py_funcname);
    return NULL;
}

static void __Pyx_AddTraceback(const char *funcname, int c_line,
                               int py_line, const char *filename) {
    PyCodeObject *py_code = 0;
    PyFrameObject *py_frame = 0;
    PyThreadState *tstate = __Pyx_PyThreadState_Current;
    PyObject *ptype, *pvalue, *ptraceback;

    if (c_line) {
        c_line = __Pyx_CLineForTraceback(tstate, c_line);
    }

    // Negate to avoid collisions between py and c lines.
    py_code = $global_code_object_cache_find(c_line ? -c_line : py_line);
    if (!py_code) {
        __Pyx_ErrFetchInState(tstate, &ptype, &pvalue, &ptraceback);
        py_code = __Pyx_CreateCodeObjectForTraceback(
            funcname, c_line, py_line, filename);
        if (!py_code) {
            /* If the code object creation fails, then we should clear the
               fetched exception references and propagate the new exception */
            Py_XDECREF(ptype);
            Py_XDECREF(pvalue);
            Py_XDECREF(ptraceback);
            goto bad;
        }
        __Pyx_ErrRestoreInState(tstate, ptype, pvalue, ptraceback);
        $global_code_object_cache_insert(c_line ? -c_line : py_line, py_code);
    }
    py_frame = PyFrame_New(
        tstate,            /*PyThreadState *tstate,*/
        py_code,           /*PyCodeObject *code,*/
        NAMED_CGLOBAL(moddict_cname),    /*PyObject *globals,*/
        0                  /*PyObject *locals*/
    );
    if (!py_frame) goto bad;
    __Pyx_PyFrame_SetLineNumber(py_frame, py_line);
    PyTraceBack_Here(py_frame);
bad:
    Py_XDECREF(py_code);
    Py_XDECREF(py_frame);
}
#endif

///////////////////////////// FloatExceptionCheck.proto ///////////////////////////

// Detect if error_value is NaN, and use a different check in that case
#define __PYX_CHECK_FLOAT_EXCEPTION(value, error_value) \
    ((error_value) == (error_value) ? \
     (value) == (error_value) : \
     (value) != (value))
