/////////////// Profile_config.proto ///////////////
//@proto_block: utility_code_proto_before_types

#ifndef CYTHON_PROFILE
#if CYTHON_COMPILING_IN_LIMITED_API || CYTHON_COMPILING_IN_PYPY
  #define CYTHON_PROFILE 0
#else
  #define CYTHON_PROFILE 1
#endif
#endif

#ifndef CYTHON_TRACE_NOGIL
  #define CYTHON_TRACE_NOGIL 0
#else
  #if CYTHON_TRACE_NOGIL && !defined(CYTHON_TRACE)
    #define CYTHON_TRACE 1
  #endif
#endif

#ifndef CYTHON_TRACE
  #define CYTHON_TRACE 0
#endif


#if CYTHON_PROFILE || CYTHON_TRACE
#if CYTHON_USE_SYS_MONITORING

    // TODO: make list of event types specific to functions.
    typedef enum {
        __Pyx_Monitoring_PY_START = 0,
        __Pyx_Monitoring_PY_RETURN,
        __Pyx_Monitoring_PY_UNWIND,
        //__Pyx_Monitoring_CALL,
        __Pyx_Monitoring_LINE,
        __Pyx_Monitoring_RAISE,
        __Pyx_Monitoring_RERAISE,
        __Pyx_Monitoring_EXCEPTION_HANDLED,
        __Pyx_Monitoring_PY_RESUME,
        __Pyx_Monitoring_PY_YIELD,
        __Pyx_Monitoring_STOP_ITERATION,
    } __Pyx_Monitoring_Event_Index;

    // Note: This is used by the generator/coroutine object struct (in its '.proto' section).
    static const unsigned char __Pyx_MonitoringEventTypes[] = {
        PY_MONITORING_EVENT_PY_START,
        PY_MONITORING_EVENT_PY_RETURN,
        PY_MONITORING_EVENT_PY_UNWIND,
        //PY_MONITORING_EVENT_CALL,
        PY_MONITORING_EVENT_LINE,
        PY_MONITORING_EVENT_RAISE,
        PY_MONITORING_EVENT_RERAISE,
        PY_MONITORING_EVENT_EXCEPTION_HANDLED,
        // generator specific:
        PY_MONITORING_EVENT_PY_RESUME,
        PY_MONITORING_EVENT_PY_YIELD,
        PY_MONITORING_EVENT_STOP_ITERATION,
    };

    #define __Pyx_MonitoringEventTypes_CyFunc_count (sizeof(__Pyx_MonitoringEventTypes) - 3)
    #define __Pyx_MonitoringEventTypes_CyGen_count (sizeof(__Pyx_MonitoringEventTypes))

#endif
#endif


/////////////// Profile.proto ///////////////
//@requires: Exceptions.c::PyErrFetchRestore
//@requires: Profile_config
//@substitute: naming

// Note that cPython ignores PyTrace_EXCEPTION,
// but maybe some other profilers don't.

#if CYTHON_TRACE
  #undef CYTHON_PROFILE_REUSE_FRAME
#endif

#if CYTHON_USE_MODULE_STATE
  #undef CYTHON_PROFILE_REUSE_CODEOBJ
  #define CYTHON_PROFILE_REUSE_CODEOBJ 0
  #undef CYTHON_PROFILE_REUSE_FRAME
#endif

#ifndef CYTHON_PROFILE_REUSE_CODEOBJ
  #define CYTHON_PROFILE_REUSE_CODEOBJ 1
#endif

#ifndef CYTHON_PROFILE_REUSE_FRAME
  #define CYTHON_PROFILE_REUSE_FRAME 0
#endif

#if CYTHON_USE_SYS_MONITORING && (CYTHON_PROFILE || CYTHON_TRACE)
  // Some types are shared across modules with an object struct layout that depends on monitoring support.
  // We therefore use a different name for the different types.
  #define __PYX_MONITORING_ABI_SUFFIX  "_mon"
#else
  #define __PYX_MONITORING_ABI_SUFFIX
#endif

#if CYTHON_PROFILE || CYTHON_TRACE

#if CYTHON_USE_SYS_MONITORING
// Use the Py3.13 monitoring C-API: https://github.com/python/cpython/issues/111997

  typedef uint64_t __pyx_monitoring_version_type;

  #define __Pyx_TraceDeclarationsFunc \
      PyObject *$frame_code_cname = NULL; \
      PyMonitoringState $monitoring_states_cname[__Pyx_MonitoringEventTypes_CyFunc_count]; \
      int __pyx_exception_already_reported = 0; \
      const int __pyx_sys_monitoring_disabled_in_parallel = 0; CYTHON_UNUSED_VAR(__pyx_sys_monitoring_disabled_in_parallel);

  #define __Pyx_TraceDeclarationsGen \
      PyObject *$frame_code_cname = Py_NewRef($generator_cname->gi_code); \
      PyMonitoringState* $monitoring_states_cname = $generator_cname->$monitoring_states_cname; \
      __pyx_monitoring_version_type $monitoring_version_cname = $generator_cname->$monitoring_version_cname; \
      int __pyx_exception_already_reported = 0; \
      const int __pyx_sys_monitoring_disabled_in_parallel = 0; CYTHON_UNUSED_VAR(__pyx_sys_monitoring_disabled_in_parallel);

  #define __Pyx_IsTracing(event_id)  ((!__pyx_sys_monitoring_disabled_in_parallel) && ($monitoring_states_cname[event_id]).active)
  #define __Pyx_TraceFrameInit(codeobj) \
      if (codeobj) $frame_code_cname = codeobj;

  #define __Pyx_TurnOffSysMonitoringInParallel \
    const int __pyx_sys_monitoring_disabled_in_parallel = 1; \
    CYTHON_UNUSED_VAR(__pyx_sys_monitoring_disabled_in_parallel);

  CYTHON_UNUSED static PyCodeObject *__Pyx_createFrameCodeObject(const char *funcname, const char *srcfile, int firstlineno); /*proto*/
  CYTHON_UNUSED static int __Pyx__TraceStartFunc(PyMonitoringState *state_array, PyObject *code_obj, int offset, int skip_event); /*proto*/
  CYTHON_UNUSED static int __Pyx__TraceStartGen(PyMonitoringState *state_array, __pyx_monitoring_version_type *version, PyObject *code_obj, int offset); /*proto*/
  CYTHON_UNUSED static int __Pyx__TraceResumeGen(PyMonitoringState *state_array, __pyx_monitoring_version_type *version, PyObject *code_obj, int offset); /*proto*/
  CYTHON_UNUSED static void __Pyx__TraceException(PyMonitoringState *monitoring_state, PyObject *code_obj, int offset, int reraised); /*proto*/

  #define __Pyx_PyMonitoring_ExitScope(nogil)              \
    if (nogil) {                                           \
        (void) __pyx_exception_already_reported;           \
        if (CYTHON_TRACE_NOGIL) {                          \
            PyGILState_STATE state = PyGILState_Ensure();  \
            PyMonitoring_ExitScope();                      \
            Py_XDECREF($frame_code_cname);                 \
            PyGILState_Release(state);                     \
        }                                                  \
    } else {                                               \
        PyMonitoring_ExitScope();                          \
        Py_XDECREF($frame_code_cname);                     \
    }


  // We check "tstate->tracing" after clearing the monitoring state to prevent re-entry while a trace function is running.
  #define __Pyx_TraceStartFunc(funcname, srcfile, firstlineno, offset, nogil, skip_event, goto_error) \
  if ((0) /* !__Pyx_IsTracing(__Pyx_Monitoring_PY_START) */); else {                         \
      int ret = 0;                                                                           \
      memset($monitoring_states_cname, 0, sizeof($monitoring_states_cname));                 \
      if (nogil) {                                                                           \
          if (CYTHON_TRACE_NOGIL) {                                                          \
              PyGILState_STATE state = PyGILState_Ensure();                                  \
              if (!__Pyx_PyThreadState_Current->tracing) {                                   \
                  if (likely($frame_code_cname)) Py_INCREF($frame_code_cname);               \
                  else $frame_code_cname = (PyObject*) __Pyx_createFrameCodeObject(funcname, srcfile, firstlineno); \
                  if (unlikely(!$frame_code_cname)) ret = -1;                                \
                  else ret = __Pyx__TraceStartFunc($monitoring_states_cname, $frame_code_cname, offset, skip_event); \
              } else $frame_code_cname = NULL;                                               \
              PyGILState_Release(state);                                                     \
          } else $frame_code_cname = NULL;                                                   \
      } else {                                                                               \
          if (!__Pyx_PyThreadState_Current->tracing) {                                       \
              if (likely($frame_code_cname)) Py_INCREF($frame_code_cname);                   \
              else $frame_code_cname = (PyObject*) __Pyx_createFrameCodeObject(funcname, srcfile, firstlineno); \
              if (unlikely(!$frame_code_cname)) ret = -1;                                    \
              else ret = __Pyx__TraceStartFunc($monitoring_states_cname, $frame_code_cname, offset, skip_event); \
          } else $frame_code_cname = NULL;                                                   \
      }                                                                                      \
      if (unlikely(ret == -1)) goto_error;                                                   \
  }

  // TODO: We should prevent tracing inside of trace functions (tstate->tracing > 0).
  // 'nogil' is obviously unused
  #define __Pyx_TraceStartGen(funcname, srcfile, firstlineno, offset, nogil, skip_event, goto_error) \
  if ((0) /* !__Pyx_IsTracing(__Pyx_Monitoring_PY_START) */); else {                         \
      int ret = __Pyx__TraceStartGen($monitoring_states_cname, &$monitoring_version_cname, $frame_code_cname, offset); \
      if (unlikely(ret == -1)) goto_error;                                                   \
  }

  // TODO: We should prevent tracing inside of trace functions (tstate->tracing > 0).
  #define __Pyx_TraceResumeGen(funcname, srcfile, firstlineno, offset, goto_error) \
  if ((0) /* !__Pyx_IsTracing(__Pyx_Monitoring_PY_RESUME) */); else {                        \
      int ret = __Pyx__TraceResumeGen($monitoring_states_cname, &$monitoring_version_cname, $frame_code_cname, offset); \
      if (unlikely(ret == -1)) goto_error;                                                   \
  }

  #define __Pyx_TraceYield(result, offset, goto_error) \
  if (!__Pyx_IsTracing(__Pyx_Monitoring_PY_YIELD)); else {                                   \
      int ret = PyMonitoring_FirePyYieldEvent(&$monitoring_states_cname[__Pyx_Monitoring_PY_RETURN], $frame_code_cname, offset, result); \
      PyMonitoring_ExitScope();                                                              \
      if (unlikely(ret == -1)) goto_error;                                                   \
  }

//   #if PY_VERSION_HEX >= 0x030d00b2
//   #define __Pyx_TraceStopIteration(value, offset, goto_error) \
//   if (!__Pyx_IsTracing(__Pyx_Monitoring_STOP_ITERATION)); else {                             \
//       int ret = PyMonitoring_FireStopIterationEvent(&$monitoring_states_cname[__Pyx_Monitoring_STOP_ITERATION], $frame_code_cname, offset, value); \
//       if (unlikely(ret == -1)) goto_error;                                                   \
//   }
//   #else
//   #define __Pyx_TraceStopIteration(value, offset, goto_error) \
//   if (!__Pyx_IsTracing(__Pyx_Monitoring_STOP_ITERATION)); else {                             \
//       PyErr_SetObject(PyExc_StopIteration, value);                                           \
//       int ret = PyMonitoring_FireStopIterationEvent(&$monitoring_states_cname[__Pyx_Monitoring_STOP_ITERATION], $frame_code_cname, offset, value); \
//       if (unlikely(ret == -1)) goto_error;                                                   \
//       PyErr_SetRaisedException(NULL);                                                        \
//   }
//   #endif

  // No error handling here since the exception path will be taken either way.
  #define __Pyx_TraceException(offset, reraised, fresh) \
  if (!__Pyx_IsTracing((reraised) ? __Pyx_Monitoring_RERAISE : __Pyx_Monitoring_RAISE)); else { \
      if (fresh || reraised || !__pyx_exception_already_reported) {                             \
          __Pyx__TraceException(&$monitoring_states_cname[(reraised) ? __Pyx_Monitoring_RERAISE : __Pyx_Monitoring_RAISE], $frame_code_cname, offset, reraised); \
      }                                                                                         \
      __pyx_exception_already_reported = 1;                                                     \
  }

  #define __Pyx_TraceExceptionDone()  __pyx_exception_already_reported = 0

  // No error handling here since the exception path we'll catch the exception next.
  #define __Pyx_TraceExceptionHandled(offset) \
  if (!__Pyx_IsTracing(__Pyx_Monitoring_EXCEPTION_HANDLED)); else {                          \
      (void) PyMonitoring_FireExceptionHandledEvent(&$monitoring_states_cname[__Pyx_Monitoring_EXCEPTION_HANDLED], $frame_code_cname, offset); \
      __pyx_exception_already_reported = 0;                                                  \
  }

  // We assume that we own a safe reference to the returned value, usually in `__pyx_r`.
  #define __Pyx_TraceReturnValue(result, offset, nogil, goto_error) \
  if (!__Pyx_IsTracing(__Pyx_Monitoring_PY_RETURN)); else {                                  \
      int ret = 0;                                                                           \
      if (nogil) {                                                                           \
          if (CYTHON_TRACE_NOGIL) {                                                          \
              PyGILState_STATE state = PyGILState_Ensure();                                  \
              ret = PyMonitoring_FirePyReturnEvent(&$monitoring_states_cname[__Pyx_Monitoring_PY_RETURN], $frame_code_cname, offset, result); \
              PyGILState_Release(state);                                                     \
          }                                                                                  \
      } else {                                                                               \
          ret = PyMonitoring_FirePyReturnEvent(&$monitoring_states_cname[__Pyx_Monitoring_PY_RETURN], $frame_code_cname, offset, result); \
      }                                                                                      \
      if (unlikely(ret == -1)) goto_error;                                                   \
  }

  #define __Pyx_TraceReturnCValue(cresult, convert_function, offset, nogil, goto_error) \
  if (!__Pyx_IsTracing(__Pyx_Monitoring_PY_RETURN)); else {                                  \
      int ret = 0;                                                                           \
      if (nogil) {                                                                           \
          if (CYTHON_TRACE_NOGIL) {                                                          \
              PyGILState_STATE state = PyGILState_Ensure();                                  \
              PyObject *pyvalue = convert_function(cresult);                                 \
              if (unlikely(!pyvalue)) {                                                      \
                  PyErr_Clear();                                                             \
                  pyvalue = Py_None; Py_INCREF(Py_None);                                     \
              }                                                                              \
              ret = PyMonitoring_FirePyReturnEvent(&$monitoring_states_cname[__Pyx_Monitoring_PY_RETURN], $frame_code_cname, offset, pyvalue); \
              Py_DECREF(pyvalue);                                                            \
              PyGILState_Release(state);                                                     \
          }                                                                                  \
      } else {                                                                               \
          PyObject *pyvalue = convert_function(cresult);                                     \
          if (unlikely(!pyvalue)) {                                                          \
              PyErr_Clear();                                                                 \
              pyvalue = Py_None; Py_INCREF(Py_None);                                         \
          }                                                                                  \
          ret = PyMonitoring_FirePyReturnEvent(&$monitoring_states_cname[__Pyx_Monitoring_PY_RETURN], $frame_code_cname, offset, pyvalue); \
          Py_DECREF(pyvalue);                                                                \
      }                                                                                      \
      if (unlikely(ret == -1)) goto_error;                                                   \
  }

  // We don't need an exception goto in __Pyx_TraceExceptionUnwind() because we are propagating an exception either way.
  #define __Pyx_TraceExceptionUnwind(offset, nogil) \
  if (!__Pyx_IsTracing(__Pyx_Monitoring_PY_UNWIND)); else {                                  \
      if (nogil) {                                                                           \
          if (CYTHON_TRACE_NOGIL) {                                                          \
              PyGILState_STATE state = PyGILState_Ensure();                                  \
              (void) PyMonitoring_FirePyUnwindEvent(&$monitoring_states_cname[__Pyx_Monitoring_PY_UNWIND], $frame_code_cname, offset); \
              PyGILState_Release(state);                                                     \
          }                                                                                  \
      } else {                                                                               \
          (void) PyMonitoring_FirePyUnwindEvent(&$monitoring_states_cname[__Pyx_Monitoring_PY_UNWIND], $frame_code_cname, offset); \
      }                                                                                      \
  }

  #if CYTHON_TRACE

  CYTHON_UNUSED static int __Pyx__TraceLine(PyMonitoringState *monitoring_state, PyObject *code_obj, int line, int offset); /*proto*/

  #define __Pyx_TraceLine(line, offset, nogil, goto_error) \
  if (!__Pyx_IsTracing(__Pyx_Monitoring_LINE)); else {                                       \
      int ret = 0;                                                                           \
      if (nogil) {                                                                           \
          if (CYTHON_TRACE_NOGIL) {                                                          \
              PyGILState_STATE state = PyGILState_Ensure();                                  \
              ret = __Pyx__TraceLine(&$monitoring_states_cname[__Pyx_Monitoring_LINE], $frame_code_cname, line, offset); \
              PyGILState_Release(state);                                                     \
          }                                                                                  \
      } else {                                                                               \
          ret = __Pyx__TraceLine(&$monitoring_states_cname[__Pyx_Monitoring_LINE], $frame_code_cname, line, offset); \
      }                                                                                      \
      if (unlikely(ret == -1)) goto_error;                                                   \
  }
  #endif

#else
// No PyMonitoring C-API (Py3.13+), use pre-Py3.12 profiling/tracing "C-API".

  #include "compile.h"
  #include "frameobject.h"
  #include "traceback.h"
#if PY_VERSION_HEX >= 0x030b00a6 && !defined(PYPY_VERSION)
  #ifndef Py_BUILD_CORE
    #define Py_BUILD_CORE 1
  #endif
  #include "internal/pycore_frame.h"
#endif

  #if CYTHON_PROFILE_REUSE_FRAME
    #define CYTHON_FRAME_MODIFIER static
    #define CYTHON_FRAME_DEL(frame)
  #else
    #define CYTHON_FRAME_MODIFIER
    #define CYTHON_FRAME_DEL(frame) Py_CLEAR(frame)
  #endif

  #if CYTHON_PROFILE_REUSE_CODEOBJ
    #define CYTHON_CODEOBJ_MODIFIER static
  #else
    #define CYTHON_CODEOBJ_MODIFIER
  #endif


  #define __Pyx_TraceDeclarationsFunc \
      CYTHON_CODEOBJ_MODIFIER PyCodeObject *$frame_code_cname = NULL;     \
      CYTHON_FRAME_MODIFIER PyFrameObject *$frame_cname = NULL;           \
      int __Pyx_use_tracing = 0;

  #define __Pyx_TraceDeclarationsGen \
      PyObject *$frame_code_cname = $generator_cname->gi_code; \
      CYTHON_FRAME_MODIFIER PyFrameObject *$frame_cname = NULL; \
      int __Pyx_use_tracing = 0;

  #define __Pyx_TraceFrameInit(codeobj)                                   \
      if (codeobj) $frame_code_cname = (PyCodeObject*) codeobj;

  #define __Pyx_PyMonitoring_ExitScope(nogil)         \
    if (!CYTHON_PROFILE_REUSE_FRAME && nogil) {       \
        PyGILState_STATE state = PyGILState_Ensure(); \
        CYTHON_FRAME_DEL($frame_cname);               \
        PyGILState_Release(state);                    \
    } else {                                          \
        CYTHON_FRAME_DEL($frame_cname);               \
    }
  #define __Pyx_TraceException(offset, reraised, fresh)  {}
  #define __Pyx_TraceExceptionHandled(offset)  {}
  #define __Pyx_TraceExceptionDone()  {}
  #define __Pyx_TurnOffSysMonitoringInParallel {} // Only needed for freethreading

#if PY_VERSION_HEX >= 0x030b00a2
  #if PY_VERSION_HEX >= 0x030C00b1
  #define __Pyx_IsTracing(tstate, check_tracing, check_funcs) \
     ((!(check_tracing) || !(tstate)->tracing) && \
         (!(check_funcs) || (tstate)->c_profilefunc || (CYTHON_TRACE && (tstate)->c_tracefunc)))
  #else
  #define __Pyx_IsTracing(tstate, check_tracing, check_funcs) \
     (unlikely((tstate)->cframe->use_tracing) && \
         (!(check_tracing) || !(tstate)->tracing) && \
         (!(check_funcs) || (tstate)->c_profilefunc || (CYTHON_TRACE && (tstate)->c_tracefunc)))
  #endif

  #define __Pyx_EnterTracing(tstate)  PyThreadState_EnterTracing(tstate)
  #define __Pyx_LeaveTracing(tstate)  PyThreadState_LeaveTracing(tstate)

#elif PY_VERSION_HEX >= 0x030a00b1
  #define __Pyx_IsTracing(tstate, check_tracing, check_funcs) \
     (unlikely((tstate)->cframe->use_tracing) && \
         (!(check_tracing) || !(tstate)->tracing) && \
         (!(check_funcs) || (tstate)->c_profilefunc || (CYTHON_TRACE && (tstate)->c_tracefunc)))

  #define __Pyx_EnterTracing(tstate) \
      do { tstate->tracing++; tstate->cframe->use_tracing = 0; } while (0)

  #define __Pyx_LeaveTracing(tstate) \
      do { \
          tstate->tracing--; \
          tstate->cframe->use_tracing = ((CYTHON_TRACE && tstate->c_tracefunc != NULL) \
                                 || tstate->c_profilefunc != NULL); \
      } while (0)

#else
  #define __Pyx_IsTracing(tstate, check_tracing, check_funcs) \
     (unlikely((tstate)->use_tracing) && \
         (!(check_tracing) || !(tstate)->tracing) && \
         (!(check_funcs) || (tstate)->c_profilefunc || (CYTHON_TRACE && (tstate)->c_tracefunc)))

  #define __Pyx_EnterTracing(tstate) \
      do { tstate->tracing++; tstate->use_tracing = 0; } while (0)

  #define __Pyx_LeaveTracing(tstate) \
      do { \
          tstate->tracing--; \
          tstate->use_tracing = ((CYTHON_TRACE && tstate->c_tracefunc != NULL) \
                                         || tstate->c_profilefunc != NULL); \
      } while (0)

#endif

  #define __Pyx_TraceStartFunc(funcname, srcfile, firstlineno, offset, nogil, skip_event, goto_error) \
  if (nogil) {                                                                           \
      if (CYTHON_TRACE_NOGIL) {                                                          \
          PyThreadState *tstate;                                                         \
          PyGILState_STATE state = PyGILState_Ensure();                                  \
          tstate = __Pyx_PyThreadState_Current;                                          \
          if (__Pyx_IsTracing(tstate, 1, 1)) {                                           \
              __Pyx_use_tracing = __Pyx_TraceSetupAndCall((PyCodeObject**)&$frame_code_cname, &$frame_cname, tstate, funcname, srcfile, firstlineno, skip_event);  \
          }                                                                              \
          PyGILState_Release(state);                                                     \
          if (unlikely(__Pyx_use_tracing < 0)) goto_error;                               \
      }                                                                                  \
  } else {                                                                               \
      PyThreadState* tstate = PyThreadState_GET();                                       \
      if (__Pyx_IsTracing(tstate, 1, 1)) {                                               \
          __Pyx_use_tracing = __Pyx_TraceSetupAndCall((PyCodeObject**)&$frame_code_cname, &$frame_cname, tstate, funcname, srcfile, firstlineno, skip_event);  \
          if (unlikely(__Pyx_use_tracing < 0)) goto_error;                               \
      }                                                                                  \
  }

  #define __Pyx_TraceStartGen __Pyx_TraceStartFunc

  #define __Pyx_TraceYield(result, offset, goto_error) \
  if (likely(!__Pyx_use_tracing)); else {                                                \
      PyThreadState* tstate = __Pyx_PyThreadState_Current;                               \
      if (__Pyx_IsTracing(tstate, 0, 0)) {                                               \
          __Pyx_call_return_trace_func(tstate, $frame_cname, (PyObject*)result);         \
      }                                                                                  \
      if ((1)); else goto_error;                                                         \
  }

  #define __Pyx_TraceResumeGen(funcname, srcfile, firstlineno, offset, goto_error) \
      __Pyx_TraceStartFunc(funcname, srcfile, firstlineno, offset, 0, 0, goto_error)

  CYTHON_UNUSED static void __Pyx_call_return_trace_func(PyThreadState *tstate, PyFrameObject *frame, PyObject *result) {
      PyObject *type, *value, *traceback;
      __Pyx_ErrFetchInState(tstate, &type, &value, &traceback);
      __Pyx_EnterTracing(tstate);
      if (CYTHON_TRACE && tstate->c_tracefunc)
          tstate->c_tracefunc(tstate->c_traceobj, frame, PyTrace_RETURN, result);
      if (tstate->c_profilefunc)
          tstate->c_profilefunc(tstate->c_profileobj, frame, PyTrace_RETURN, result);
      __Pyx_LeaveTracing(tstate);
      __Pyx_ErrRestoreInState(tstate, type, value, traceback);
  }

  #define __Pyx_TraceReturnValue(result, offset, nogil, goto_error) \
  if (likely(!__Pyx_use_tracing)); else {                                                 \
      if (nogil) {                                                                        \
          if (CYTHON_TRACE_NOGIL) {                                                       \
              PyThreadState *tstate;                                                      \
              PyGILState_STATE state = PyGILState_Ensure();                               \
              tstate = __Pyx_PyThreadState_Current;                                       \
              if (__Pyx_IsTracing(tstate, 0, 0)) {                                        \
                  __Pyx_call_return_trace_func(tstate, $frame_cname, (PyObject*)result);  \
              }                                                                           \
              PyGILState_Release(state);                                                  \
          }                                                                               \
      } else {                                                                            \
          PyThreadState* tstate = __Pyx_PyThreadState_Current;                            \
          if (__Pyx_IsTracing(tstate, 0, 0)) {                                            \
              __Pyx_call_return_trace_func(tstate, $frame_cname, (PyObject*)result);      \
          }                                                                               \
      }                                                                                   \
      if ((1)); else goto_error;                                                          \
  }

  #define __Pyx_TraceReturnCValue(cresult, convert_function, offset, nogil, goto_error) \
  if (likely(!__Pyx_use_tracing)); else {                                                 \
      if (nogil) {                                                                        \
          if (CYTHON_TRACE_NOGIL) {                                                       \
              PyThreadState *tstate;                                                      \
              PyGILState_STATE state = PyGILState_Ensure();                               \
              tstate = __Pyx_PyThreadState_Current;                                       \
              if (__Pyx_IsTracing(tstate, 0, 0)) {                                        \
                  PyObject *pyvalue = convert_function(cresult);                          \
                  if (unlikely(!pyvalue)) {                                               \
                    PyErr_Clear();                                                        \
                    pyvalue = Py_None; Py_INCREF(Py_None);                                \
                  }                                                                       \
                  __Pyx_call_return_trace_func(tstate, $frame_cname, pyvalue);            \
                  Py_DECREF(pyvalue);                                                     \
              }                                                                           \
              PyGILState_Release(state);                                                  \
          }                                                                               \
      } else {                                                                            \
          PyThreadState* tstate = __Pyx_PyThreadState_Current;                            \
          if (__Pyx_IsTracing(tstate, 0, 0)) {                                            \
              PyObject *pyvalue = convert_function(cresult);                              \
              if (unlikely(!pyvalue)) {                                                   \
                  PyErr_Clear();                                                          \
                  pyvalue = Py_None; Py_INCREF(Py_None);                                  \
              }                                                                           \
              __Pyx_call_return_trace_func(tstate, $frame_cname, pyvalue);                \
              Py_DECREF(pyvalue);                                                         \
          }                                                                               \
      }                                                                                   \
      if ((1)); else goto_error;                                                          \
  }

  #define __Pyx_TraceExceptionUnwind(offset, nogil) \
  if (likely(!__Pyx_use_tracing)); else {                                                 \
      if (nogil) {                                                                        \
          if (CYTHON_TRACE_NOGIL) {                                                       \
              PyThreadState *tstate;                                                      \
              PyGILState_STATE state = PyGILState_Ensure();                               \
              tstate = __Pyx_PyThreadState_Current;                                       \
              if (__Pyx_IsTracing(tstate, 0, 0)) {                                        \
                  __Pyx_call_return_trace_func(tstate, $frame_cname, Py_None);            \
              }                                                                           \
              PyGILState_Release(state);                                                  \
          }                                                                               \
      } else {                                                                            \
          PyThreadState* tstate = __Pyx_PyThreadState_Current;                            \
          if (__Pyx_IsTracing(tstate, 0, 0)) {                                            \
              __Pyx_call_return_trace_func(tstate, $frame_cname, Py_None);                \
          }                                                                               \
      }                                                                                   \
  }

  static int __Pyx_TraceSetupAndCall(PyCodeObject** code, PyFrameObject** frame, PyThreadState* tstate, const char *funcname, const char *srcfile, int firstlineno, int skip_event); /*proto*/

#if CYTHON_TRACE
  CYTHON_UNUSED static int __Pyx_call_line_trace_func(PyThreadState *tstate, PyFrameObject *frame, int line);/*proto*/

  #define __Pyx_TraceLine(line, offset, nogil, goto_error)                                 \
  if (likely(!__Pyx_use_tracing)); else {                                                  \
      int ret = 0;                                                                         \
      if (nogil) {                                                                         \
          if (CYTHON_TRACE_NOGIL) {                                                        \
              PyThreadState *tstate;                                                       \
              PyGILState_STATE state = __Pyx_PyGILState_Ensure();                          \
              tstate = __Pyx_PyThreadState_Current;                                        \
              if (__Pyx_IsTracing(tstate, 0, 0) && tstate->c_tracefunc && $frame_cname->f_trace) { \
                  ret = __Pyx_call_line_trace_func(tstate, $frame_cname, line);            \
              }                                                                            \
              __Pyx_PyGILState_Release(state);                                             \
          }                                                                                \
      } else {                                                                             \
          PyThreadState* tstate = __Pyx_PyThreadState_Current;                             \
          if (__Pyx_IsTracing(tstate, 0, 0) && tstate->c_tracefunc && $frame_cname->f_trace) { \
              ret = __Pyx_call_line_trace_func(tstate, $frame_cname, line);                \
          }                                                                                \
      }                                                                                    \
      if (unlikely(ret)) goto_error;                                                       \
  }
#endif

// End of pre-monitoring implementation (Py<3.12)
#endif

#else
// ! (CYTHON_TRACE || CYTHON_PROFILE)

  #define __Pyx_TraceDeclarationsFunc
  #define __Pyx_TraceDeclarationsGen
  #define __Pyx_TraceExceptionDone()  {}
  #define __Pyx_TraceFrameInit(codeobj)  {}
  #define __Pyx_TurnOffSysMonitoringInParallel {}
  #define __Pyx_PyMonitoring_ExitScope(nogil)  {}
  #define __Pyx_TraceException(offset, reraised, fresh)  {}
  #define __Pyx_TraceExceptionUnwind(offset, nogil)  {}
  #define __Pyx_TraceExceptionHandled(offset)  {}
  // mark error label as used to avoid compiler warnings
  #define __Pyx_TraceStartFunc(funcname, srcfile, firstlineno, offset, nogil, skip_event, goto_error)   if ((1)); else goto_error;
  #define __Pyx_TraceStartGen __Pyx_TraceStartFunc
  #define __Pyx_TraceResumeGen(funcname, srcfile, firstlineno, offset, goto_error)   if ((1)); else goto_error;
  #define __Pyx_TraceYield(result, offset, goto_error)   if ((1)); else goto_error;
  #define __Pyx_TraceReturnValue(result, offset, nogil, goto_error) \
      if ((1)); else goto_error;
  #define __Pyx_TraceReturnCValue(cresult, convert_function, offset, nogil, goto_error) \
      if ((1)); else { (void) convert_function; goto_error }

#endif /* CYTHON_PROFILE || CYTHON_TRACE */

#if !CYTHON_TRACE
  // mark error label as used to avoid compiler warnings
  #define __Pyx_TraceLine(line, offset, nogil, goto_error)   if ((1)); else goto_error;
#endif

/////////////// Profile ///////////////

#if CYTHON_PROFILE || CYTHON_TRACE

#if CYTHON_TRACE && !CYTHON_USE_SYS_MONITORING
static int __Pyx_call_line_trace_func(PyThreadState *tstate, PyFrameObject *frame, int line) {
    // see call_trace_protected() in CPython's ceval.c
    int ret;
    PyObject *type, *value, *traceback;
    __Pyx_ErrFetchInState(tstate, &type, &value, &traceback);
    __Pyx_PyFrame_SetLineNumber(frame, line);
    __Pyx_EnterTracing(tstate);

    ret = tstate->c_tracefunc(tstate->c_traceobj, frame, PyTrace_LINE, NULL);

    __Pyx_LeaveTracing(tstate);
    if (likely(!ret)) {
        __Pyx_ErrRestoreInState(tstate, type, value, traceback);
    } else {
        Py_XDECREF(type);
        Py_XDECREF(value);
        Py_XDECREF(traceback);
    }
    return ret;
}
#endif

CYTHON_UNUSED static PyCodeObject *__Pyx_createFrameCodeObject(const char *funcname, const char *srcfile, int firstlineno) {
    PyCodeObject *py_code = PyCode_NewEmpty(srcfile, funcname, firstlineno);
    // make CPython use a fresh dict for "f_locals" at need (see GH #1836)
    if (likely(py_code)) {
        py_code->co_flags |= CO_OPTIMIZED | CO_NEWLOCALS;
    }
    return py_code;
}

#if CYTHON_USE_SYS_MONITORING

CYTHON_UNUSED static int __Pyx__TraceStartFunc(PyMonitoringState *state_array, PyObject *code_obj, int offset, int skip_event) {
    int ret;
    __pyx_monitoring_version_type version = 0;
    ret = PyMonitoring_EnterScope(state_array, &version, __Pyx_MonitoringEventTypes, __Pyx_MonitoringEventTypes_CyFunc_count);
    if (unlikely(ret == -1)) return -1;
    return skip_event ? 0 : PyMonitoring_FirePyStartEvent(&state_array[__Pyx_Monitoring_PY_START], code_obj, offset);
}

CYTHON_UNUSED static int __Pyx__TraceStartGen(PyMonitoringState *state_array, __pyx_monitoring_version_type *version, PyObject *code_obj, int offset) {
    int ret;
    ret = PyMonitoring_EnterScope(state_array, version, __Pyx_MonitoringEventTypes, __Pyx_MonitoringEventTypes_CyGen_count);
    if (unlikely(ret == -1)) return -1;
    return PyMonitoring_FirePyStartEvent(&state_array[__Pyx_Monitoring_PY_START], code_obj, offset);
}

CYTHON_UNUSED static int __Pyx__TraceResumeGen(PyMonitoringState *state_array, __pyx_monitoring_version_type *version, PyObject *code_obj, int offset) {
    int ret;
    ret = PyMonitoring_EnterScope(state_array, version, __Pyx_MonitoringEventTypes, __Pyx_MonitoringEventTypes_CyGen_count);
    if (unlikely(ret == -1)) return -1;
    return PyMonitoring_FirePyResumeEvent(&state_array[__Pyx_Monitoring_PY_RESUME], code_obj, offset);
}

CYTHON_UNUSED static void __Pyx__TraceException(PyMonitoringState *monitoring_state, PyObject *code_obj, int offset, int reraised) {
    if (reraised) {
        (void) PyMonitoring_FireReraiseEvent(monitoring_state, code_obj, offset);
    } else {
        (void) PyMonitoring_FireRaiseEvent(monitoring_state, code_obj, offset);
    }
}

#if CYTHON_TRACE
CYTHON_UNUSED static int __Pyx__TraceLine(PyMonitoringState *monitoring_state, PyObject *code_obj, int line, int offset) {
    int ret;
    PyObject *exc = PyErr_GetRaisedException();
    ret = PyMonitoring_FireLineEvent(monitoring_state, code_obj, offset, line);
    if (exc) PyErr_SetRaisedException(exc);
    return ret;
}
#endif

#else
// pre sys.monitoring

static int __Pyx_TraceSetupAndCall(PyCodeObject** code,
                                   PyFrameObject** frame,
                                   PyThreadState* tstate,
                                   const char *funcname,
                                   const char *srcfile,
                                   int firstlineno,
                                   int skip_event) {
    if (*frame == NULL || !CYTHON_PROFILE_REUSE_FRAME) {
        int needs_new_code_obj = (*code == NULL);
        if (needs_new_code_obj) {
            *code = __Pyx_createFrameCodeObject(funcname, srcfile, firstlineno);
            if (*code == NULL) return 0;
        }
        *frame = PyFrame_New(
            tstate,                          /*PyThreadState *tstate*/
            *code,                           /*PyCodeObject *code*/
            NAMED_CGLOBAL(moddict_cname),    /*PyObject *globals*/
            0                                /*PyObject *locals*/
        );
        if (needs_new_code_obj && !CYTHON_PROFILE_REUSE_CODEOBJ)
            Py_CLEAR(*code); // otherwise the reference is owned externally
        if (*frame == NULL) return 0;
        if (CYTHON_TRACE && (*frame)->f_trace == NULL) {
            // this enables "f_lineno" lookup, at least in CPython ...
            Py_INCREF(Py_None);
            (*frame)->f_trace = Py_None;
        }
    }

    if (!skip_event) {
        PyObject *type, *value, *traceback;
        int retval = 1;
        __Pyx_PyFrame_SetLineNumber(*frame, firstlineno);

        __Pyx_EnterTracing(tstate);
        __Pyx_ErrFetchInState(tstate, &type, &value, &traceback);

        #if CYTHON_TRACE
        if (tstate->c_tracefunc)
            retval = tstate->c_tracefunc(tstate->c_traceobj, *frame, PyTrace_CALL, NULL) == 0;
        if (retval && tstate->c_profilefunc)
        #endif
            retval = tstate->c_profilefunc(tstate->c_profileobj, *frame, PyTrace_CALL, NULL) == 0;

        __Pyx_LeaveTracing(tstate);
        if (unlikely(!retval)) {
            Py_XDECREF(type);
            Py_XDECREF(value);
            Py_XDECREF(traceback);
            return -1;
        }

        __Pyx_ErrRestoreInState(tstate, type, value, traceback);
    }

    return __Pyx_IsTracing(tstate, 0, 0);
}

#endif
#endif /* CYTHON_PROFILE */
