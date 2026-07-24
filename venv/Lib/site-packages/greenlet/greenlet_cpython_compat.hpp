/* -*- indent-tabs-mode: nil; tab-width: 4; -*- */
#ifndef GREENLET_CPYTHON_COMPAT_H
#define GREENLET_CPYTHON_COMPAT_H

/**
 * Helpers for compatibility with multiple versions of CPython.
 */

#define PY_SSIZE_T_CLEAN
#include "Python.h"

/*
Python 3.10 beta 1 changed tstate->use_tracing to a nested cframe member.
See https://github.com/python/cpython/pull/25276
We have to save and restore this as well.

Python 3.13 removed PyThreadState.cframe (GH-108035).
*/
#if PY_VERSION_HEX < 0x30D0000
#    define GREENLET_USE_CFRAME 1
#else
#    define GREENLET_USE_CFRAME 0
#endif


#if PY_VERSION_HEX >= 0x30B00A4
/*
Greenlet won't compile on anything older than Python 3.11 alpha 4 (see
https://bugs.python.org/issue46090). Summary of breaking internal changes:
- Python 3.11 alpha 1 changed how frame objects are represented internally.
  - https://github.com/python/cpython/pull/30122
- Python 3.11 alpha 3 changed how recursion limits are stored.
  - https://github.com/python/cpython/pull/29524
- Python 3.11 alpha 4 changed how exception state is stored. It also includes a
  change to help greenlet save and restore the interpreter frame "data stack".
  - https://github.com/python/cpython/pull/30122
  - https://github.com/python/cpython/pull/30234
*/
#    define GREENLET_PY311 1
#else
#    define GREENLET_PY311 0
#endif


#if PY_VERSION_HEX >= 0x30C0000
#    define GREENLET_PY312 1
#else
#    define GREENLET_PY312 0
#endif

#if PY_VERSION_HEX >= 0x30D0000
#    define GREENLET_PY313 1
#else
#    define GREENLET_PY313 0
#endif

#if PY_VERSION_HEX >= 0x30E0000
#    define GREENLET_PY314 1
#else
#    define GREENLET_PY314 0
#endif

#if PY_VERSION_HEX >= 0x30F0000
#    define GREENLET_PY315 1
#else
#    define GREENLET_PY315 0
#endif


/* _Py_DEC_REFTOTAL macro has been removed from Python 3.9 by:
  https://github.com/python/cpython/commit/49932fec62c616ec88da52642339d83ae719e924

  The symbol we use to replace it was removed by at least 3.12.
*/
#if defined(Py_REF_DEBUG) && !GREENLET_PY312
#    define GREENLET_Py_DEC_REFTOTAL _Py_RefTotal--
#else
#    define GREENLET_Py_DEC_REFTOTAL
#endif


#define G_TPFLAGS_DEFAULT Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_VERSION_TAG | Py_TPFLAGS_HAVE_GC


// bpo-43760 added PyThreadState_EnterTracing() to Python 3.11.0a2
#if PY_VERSION_HEX < 0x030B00A2 && !defined(PYPY_VERSION)
static inline void PyThreadState_EnterTracing(PyThreadState *tstate)
{
    tstate->tracing++;
    tstate->cframe->use_tracing = 0;
}
#endif

// bpo-43760 added PyThreadState_LeaveTracing() to Python 3.11.0a2
#if PY_VERSION_HEX < 0x030B00A2 && !defined(PYPY_VERSION)
static inline void PyThreadState_LeaveTracing(PyThreadState *tstate)
{
    tstate->tracing--;
    int use_tracing = (tstate->c_tracefunc != NULL
                       || tstate->c_profilefunc != NULL);
    tstate->cframe->use_tracing = use_tracing;
}
#endif

#if !defined(Py_C_RECURSION_LIMIT) && defined(C_RECURSION_LIMIT)
#  define Py_C_RECURSION_LIMIT C_RECURSION_LIMIT
#endif

// Py_IsFinalizing() became a public API in Python 3.13. Map it to the
// private _Py_IsFinalizing() on older versions so all call sites can
// use the standard name. Remove this once greenlet drops support for
// Python < 3.13.
#if !GREENLET_PY313
#  define Py_IsFinalizing() _Py_IsFinalizing()
#endif

#endif /* GREENLET_CPYTHON_COMPAT_H */
