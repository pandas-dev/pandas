#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <time.h>
#include <stdint.h>
#include "librt_time.h"
#include "pythoncapi_compat.h"
#include "mypyc_util.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

// Internal function that returns a C double for mypyc primitives
// Returns high-precision time in seconds (like time.time())
static double
time_time_internal(void) {
#ifdef _WIN32
    // Windows: Use GetSystemTimePreciseAsFileTime for ~100ns precision
    FILETIME ft;
    ULARGE_INTEGER large;

    GetSystemTimePreciseAsFileTime(&ft);
    large.LowPart = ft.dwLowDateTime;
    large.HighPart = ft.dwHighDateTime;

    // Windows FILETIME is 100-nanosecond intervals since January 1, 1601
    // 116444736000000000 = number of 100-ns intervals between 1601 and 1970
    // Convert directly to seconds: 100ns * 1e-9 = 1e-7
    int64_t intervals = large.QuadPart - 116444736000000000LL;
    return (double)intervals * 1e-7;

#else  // Unix-like systems (Linux, macOS, BSD, etc.)

    // Try clock_gettime(CLOCK_REALTIME) for nanosecond precision
    // This is available on POSIX.1-2001 and later (widely available on modern systems)
#if defined(_POSIX_TIMERS) && _POSIX_TIMERS > 0
    struct timespec ts;
    if (clock_gettime(CLOCK_REALTIME, &ts) == 0) {
        // Convert seconds and nanoseconds separately to avoid large integer operations
        return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
    }
    // Fall through to gettimeofday if clock_gettime failed
#endif

    // Fallback: gettimeofday for microsecond precision
    // This is widely available (POSIX.1-2001, BSD, etc.)
    struct timeval tv;
    if (unlikely(gettimeofday(&tv, NULL) != 0)) {
        PyErr_SetFromErrno(PyExc_OSError);
        return CPY_FLOAT_ERROR;
    }

    // Convert seconds and microseconds separately to avoid large integer operations
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
#endif
}

// Wrapper function for normal Python extension usage
static PyObject*
time_time(PyObject *self, PyObject *const *args, size_t nargs) {
    if (nargs != 0) {
        PyErr_SetString(PyExc_TypeError, "time() takes no arguments");
        return NULL;
    }

    double result = time_time_internal();
    if (result == CPY_FLOAT_ERROR) {
        return NULL;
    }
    return PyFloat_FromDouble(result);
}

static PyMethodDef librt_time_module_methods[] = {
    {"time", (PyCFunction)time_time, METH_FASTCALL,
     PyDoc_STR("Return the current time in seconds since the Unix epoch as a floating point number.")},
    {NULL, NULL, 0, NULL}
};

static int
time_abi_version(void) {
    return LIBRT_TIME_ABI_VERSION;
}

static int
time_api_version(void) {
    return LIBRT_TIME_API_VERSION;
}

static int
librt_time_module_exec(PyObject *m)
{
    // Export mypyc internal C API via capsule
    static void *time_api[LIBRT_TIME_API_LEN] = {
        (void *)time_abi_version,
        (void *)time_api_version,
        (void *)time_time_internal,
    };
    PyObject *c_api_object = PyCapsule_New((void *)time_api, "librt.time._C_API", NULL);
    if (PyModule_Add(m, "_C_API", c_api_object) < 0) {
        return -1;
    }
    return 0;
}

static PyModuleDef_Slot librt_time_module_slots[] = {
    {Py_mod_exec, librt_time_module_exec},
#ifdef Py_MOD_GIL_NOT_USED
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};

static PyModuleDef librt_time_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "time",
    .m_doc = "Fast time() function optimized for mypyc",
    .m_size = 0,
    .m_methods = librt_time_module_methods,
    .m_slots = librt_time_module_slots,
};

PyMODINIT_FUNC
PyInit_time(void)
{
    return PyModuleDef_Init(&librt_time_module);
}
