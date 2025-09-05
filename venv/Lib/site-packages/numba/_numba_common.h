#ifndef NUMBA_COMMON_H_
#define NUMBA_COMMON_H_

/* __has_attribute() is a clang / gcc-5 macro */
#ifndef __has_attribute
#   define __has_attribute(x) 0
#endif

/* This attribute marks symbols that can be shared across C objects
 * but are not exposed outside of a shared library or executable.
 * Note this is default behaviour for global symbols under Windows.
 */
#if defined(_MSC_VER)
    #define VISIBILITY_HIDDEN
    #define VISIBILITY_GLOBAL __declspec(dllexport)
#elif (__has_attribute(visibility) || (defined(__GNUC__) && __GNUC__ >= 4))
    #define VISIBILITY_HIDDEN __attribute__ ((visibility("hidden")))
    #define VISIBILITY_GLOBAL __attribute__ ((visibility("default")))
#else
    #define VISIBILITY_HIDDEN
    #define VISIBILITY_GLOBAL
#endif

/*
 * Numba's version of the PyArray_DescrCheck macro from NumPy, use it as a
 * direct replacement of NumPy's PyArray_DescrCheck to ensure binary
 * compatibility.
 *
 * Details of why this is needed:
 * NumPy 1.18 changed the definition of the PyArray_DescrCheck macro here:
 * https://github.com/numpy/numpy/commit/6108b5d1e138d07e3c9f2a4e3b1933749ad0e698
 * the result of this being that building against NumPy <1.18 would prevent
 * Numba running against NumPy >= 1.20 as noted here:
 * https://github.com/numba/numba/issues/6041#issuecomment-665132199
 *
 * This macro definition is copied from:
 * https://github.com/numpy/numpy/commit/6108b5d1e138d07e3c9f2a4e3b1933749ad0e698#diff-ad2213da23136c5fc5883d9eb2d88666R26
 *
 * NOTE: This is the NumPy 1.18 and above version of the macro.
 */
#define NUMBA_PyArray_DescrCheck(op) PyObject_TypeCheck(op, &PyArrayDescr_Type)

#endif /* NUMBA_COMMON_H_ */
