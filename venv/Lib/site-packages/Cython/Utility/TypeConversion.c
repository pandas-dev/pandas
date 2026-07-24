/////////////// TypeConversions.proto ///////////////
//@requires: StringTools.c::IncludeStringH

/* Type Conversion Predeclarations */

#define __Pyx_uchar_cast(c) ((unsigned char)c)
#define __Pyx_long_cast(x) ((long)x)

#define __Pyx_fits_Py_ssize_t(v, type, is_signed)  (    \
    (sizeof(type) < sizeof(Py_ssize_t))  ||             \
    (sizeof(type) > sizeof(Py_ssize_t) &&               \
          likely(v < (type)PY_SSIZE_T_MAX ||            \
                 v == (type)PY_SSIZE_T_MAX)  &&         \
          (!is_signed || likely(v > (type)PY_SSIZE_T_MIN ||       \
                                v == (type)PY_SSIZE_T_MIN)))  ||  \
    (sizeof(type) == sizeof(Py_ssize_t) &&              \
          (is_signed || likely(v < (type)PY_SSIZE_T_MAX ||        \
                               v == (type)PY_SSIZE_T_MAX)))  )

static CYTHON_INLINE int __Pyx_is_valid_index(Py_ssize_t i, Py_ssize_t limit) {
    // Optimisation from Section 14.2 "Bounds Checking" in
    //   https://www.agner.org/optimize/optimizing_cpp.pdf
    // See https://bugs.python.org/issue28397
    // The cast to unsigned effectively tests for "0 <= i < limit".
    return (size_t) i < (size_t) limit;
}

// fast and unsafe abs(Py_ssize_t) that ignores the overflow for (-PY_SSIZE_T_MAX-1)
#if defined (__cplusplus) && __cplusplus >= 201103L
    #include <cstdlib>
    #define __Pyx_sst_abs(value) std::abs(value)
#elif SIZEOF_INT >= SIZEOF_SIZE_T
    #define __Pyx_sst_abs(value) abs(value)
#elif SIZEOF_LONG >= SIZEOF_SIZE_T
    #define __Pyx_sst_abs(value) labs(value)
#elif defined (_MSC_VER)
    // abs() is defined for long, but 64-bits type on MSVC is long long.
    // Use MS-specific _abs64 instead.
    #define __Pyx_sst_abs(value) ((Py_ssize_t)_abs64(value))
#elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #define __Pyx_sst_abs(value) llabs(value)
#elif defined (__GNUC__)
    // gcc or clang on 64 bit windows.
    #define __Pyx_sst_abs(value) __builtin_llabs(value)
#else
    #define __Pyx_sst_abs(value) ((value<0) ? -value : value)
#endif

static CYTHON_INLINE Py_ssize_t __Pyx_ssize_strlen(const char *s);
static CYTHON_INLINE const char* __Pyx_PyObject_AsString(PyObject*);
static CYTHON_INLINE const char* __Pyx_PyObject_AsStringAndSize(PyObject*, Py_ssize_t* length);

static CYTHON_INLINE PyObject* __Pyx_PyByteArray_FromString(const char*);
#define __Pyx_PyByteArray_FromStringAndSize(s, l) PyByteArray_FromStringAndSize((const char*)s, l)
#define __Pyx_PyBytes_FromString        PyBytes_FromString
#define __Pyx_PyBytes_FromStringAndSize PyBytes_FromStringAndSize
static CYTHON_INLINE PyObject* __Pyx_PyUnicode_FromString(const char*);

#if CYTHON_ASSUME_SAFE_MACROS
    #define __Pyx_PyBytes_AsWritableString(s)     ((char*) PyBytes_AS_STRING(s))
    #define __Pyx_PyBytes_AsWritableSString(s)    ((signed char*) PyBytes_AS_STRING(s))
    #define __Pyx_PyBytes_AsWritableUString(s)    ((unsigned char*) PyBytes_AS_STRING(s))
    #define __Pyx_PyBytes_AsString(s)     ((const char*) PyBytes_AS_STRING(s))
    #define __Pyx_PyBytes_AsSString(s)    ((const signed char*) PyBytes_AS_STRING(s))
    #define __Pyx_PyBytes_AsUString(s)    ((const unsigned char*) PyBytes_AS_STRING(s))
    #define __Pyx_PyByteArray_AsString(s) PyByteArray_AS_STRING(s)
#else
    #define __Pyx_PyBytes_AsWritableString(s)     ((char*) PyBytes_AsString(s))
    #define __Pyx_PyBytes_AsWritableSString(s)    ((signed char*) PyBytes_AsString(s))
    #define __Pyx_PyBytes_AsWritableUString(s)    ((unsigned char*) PyBytes_AsString(s))
    #define __Pyx_PyBytes_AsString(s)     ((const char*) PyBytes_AsString(s))
    #define __Pyx_PyBytes_AsSString(s)    ((const signed char*) PyBytes_AsString(s))
    #define __Pyx_PyBytes_AsUString(s)    ((const unsigned char*) PyBytes_AsString(s))
    #define __Pyx_PyByteArray_AsString(s) PyByteArray_AsString(s)
#endif
#define __Pyx_PyObject_AsWritableString(s)    ((char*)(__pyx_uintptr_t) __Pyx_PyObject_AsString(s))
#define __Pyx_PyObject_AsWritableSString(s)    ((signed char*)(__pyx_uintptr_t) __Pyx_PyObject_AsString(s))
#define __Pyx_PyObject_AsWritableUString(s)    ((unsigned char*)(__pyx_uintptr_t) __Pyx_PyObject_AsString(s))
#define __Pyx_PyObject_AsSString(s)    ((const signed char*) __Pyx_PyObject_AsString(s))
#define __Pyx_PyObject_AsUString(s)    ((const unsigned char*) __Pyx_PyObject_AsString(s))
#define __Pyx_PyObject_FromCString(s)  __Pyx_PyObject_FromString((const char*)s)
#define __Pyx_PyBytes_FromCString(s)   __Pyx_PyBytes_FromString((const char*)s)
#define __Pyx_PyByteArray_FromCString(s)   __Pyx_PyByteArray_FromString((const char*)s)
#define __Pyx_PyUnicode_FromCString(s) __Pyx_PyUnicode_FromString((const char*)s)
#define __Pyx_PyUnicode_FromOrdinal(o)       PyUnicode_FromOrdinal((int)o)
#define __Pyx_PyUnicode_AsUnicode            PyUnicode_AsUnicode

static CYTHON_INLINE PyObject *__Pyx_NewRef(PyObject *obj) {
#if CYTHON_COMPILING_IN_CPYTHON && PY_VERSION_HEX >= 0x030a0000 || defined(Py_NewRef)
    return Py_NewRef(obj);
#else
    Py_INCREF(obj);
    return obj;
#endif
}

static CYTHON_INLINE PyObject *__Pyx_XNewRef(PyObject *obj) {
#if CYTHON_COMPILING_IN_CPYTHON && PY_VERSION_HEX >= 0x030a0000 || defined(Py_XNewRef)
    return Py_XNewRef(obj);
#else
    Py_XINCREF(obj);
    return obj;
#endif
}

static CYTHON_INLINE PyObject *__Pyx_Owned_Py_None(int b);
static CYTHON_INLINE PyObject * __Pyx_PyBool_FromLong(long b);
static CYTHON_INLINE int __Pyx_PyObject_IsTrue(PyObject*);
static CYTHON_INLINE int __Pyx_PyObject_IsTrueAndDecref(PyObject*);
static CYTHON_INLINE PyObject* __Pyx_PyNumber_Long(PyObject* x);

#define __Pyx_PySequence_Tuple(obj) \
    (likely(PyTuple_CheckExact(obj)) ? __Pyx_NewRef(obj) : PySequence_Tuple(obj))

static CYTHON_INLINE Py_ssize_t __Pyx_PyIndex_AsSsize_t(PyObject*);
static CYTHON_INLINE PyObject * __Pyx_PyLong_FromSize_t(size_t);
static CYTHON_INLINE Py_hash_t __Pyx_PyIndex_AsHash_t(PyObject*);

#if CYTHON_ASSUME_SAFE_MACROS
#define __Pyx_PyFloat_AsDouble(x) (PyFloat_CheckExact(x) ? PyFloat_AS_DOUBLE(x) : PyFloat_AsDouble(x))
#define __Pyx_PyFloat_AS_DOUBLE(x) PyFloat_AS_DOUBLE(x)
#else
#define __Pyx_PyFloat_AsDouble(x) PyFloat_AsDouble(x)
#define __Pyx_PyFloat_AS_DOUBLE(x) PyFloat_AsDouble(x)
#endif
#define __Pyx_PyFloat_AsFloat(x) ((float) __Pyx_PyFloat_AsDouble(x))

// We call this __Pyx_PyNumber_Int() since it's the equivalent of the int() function call.
#define __Pyx_PyNumber_Int(x) (PyLong_CheckExact(x) ? __Pyx_NewRef(x) : PyNumber_Long(x))
// __Pyx_PyNumber_Float is now in its own section since it has dependencies (needed to make
// string conversion work the same in all circumstances).

#if CYTHON_USE_PYLONG_INTERNALS
  #if PY_VERSION_HEX >= 0x030C00A7
  #ifndef _PyLong_SIGN_MASK
    #define _PyLong_SIGN_MASK 3
  #endif
  #ifndef _PyLong_NON_SIZE_BITS
    #define _PyLong_NON_SIZE_BITS 3
  #endif
  #define __Pyx_PyLong_Sign(x)  (((PyLongObject*)x)->long_value.lv_tag & _PyLong_SIGN_MASK)
  #define __Pyx_PyLong_IsNeg(x)  ((__Pyx_PyLong_Sign(x) & 2) != 0)
  #define __Pyx_PyLong_IsNonNeg(x)  (!__Pyx_PyLong_IsNeg(x))
  #define __Pyx_PyLong_IsZero(x)  (__Pyx_PyLong_Sign(x) & 1)
  #define __Pyx_PyLong_IsPos(x)  (__Pyx_PyLong_Sign(x) == 0)
  #define __Pyx_PyLong_CompactValueUnsigned(x)  (__Pyx_PyLong_Digits(x)[0])
  #define __Pyx_PyLong_DigitCount(x)  ((Py_ssize_t) (((PyLongObject*)x)->long_value.lv_tag >> _PyLong_NON_SIZE_BITS))
  #define __Pyx_PyLong_SignedDigitCount(x)  \
        ((1 - (Py_ssize_t) __Pyx_PyLong_Sign(x)) * __Pyx_PyLong_DigitCount(x))

  #if defined(PyUnstable_Long_IsCompact) && defined(PyUnstable_Long_CompactValue)
    #define __Pyx_PyLong_IsCompact(x)     PyUnstable_Long_IsCompact((PyLongObject*) x)
    #define __Pyx_PyLong_CompactValue(x)  PyUnstable_Long_CompactValue((PyLongObject*) x)
  #else
    #define __Pyx_PyLong_IsCompact(x)     (((PyLongObject*)x)->long_value.lv_tag < (2 << _PyLong_NON_SIZE_BITS))
    #define __Pyx_PyLong_CompactValue(x)  ((1 - (Py_ssize_t) __Pyx_PyLong_Sign(x)) * (Py_ssize_t) __Pyx_PyLong_Digits(x)[0])
  #endif

  // CPython 3.12 requires C99, which defines 'size_t' (but not 'ssize_t')
  typedef Py_ssize_t  __Pyx_compact_pylong;
  typedef size_t  __Pyx_compact_upylong;

  #else  /* Py < 3.12 */
  #define __Pyx_PyLong_IsNeg(x)  (Py_SIZE(x) < 0)
  #define __Pyx_PyLong_IsNonNeg(x)  (Py_SIZE(x) >= 0)
  #define __Pyx_PyLong_IsZero(x)  (Py_SIZE(x) == 0)
  #define __Pyx_PyLong_IsPos(x)  (Py_SIZE(x) > 0)
  #define __Pyx_PyLong_CompactValueUnsigned(x)  ((Py_SIZE(x) == 0) ? 0 : __Pyx_PyLong_Digits(x)[0])
  #define __Pyx_PyLong_DigitCount(x)  __Pyx_sst_abs(Py_SIZE(x))
  #define __Pyx_PyLong_SignedDigitCount(x)  Py_SIZE(x)

  #define __Pyx_PyLong_IsCompact(x)  (Py_SIZE(x) == 0 || Py_SIZE(x) == 1 || Py_SIZE(x) == -1)
  #define __Pyx_PyLong_CompactValue(x)  \
        ((Py_SIZE(x) == 0) ? (sdigit) 0 : ((Py_SIZE(x) < 0) ? -(sdigit)__Pyx_PyLong_Digits(x)[0] : (sdigit)__Pyx_PyLong_Digits(x)[0]))

  typedef sdigit  __Pyx_compact_pylong;
  typedef digit  __Pyx_compact_upylong;
  #endif

  #if PY_VERSION_HEX >= 0x030C00A5
  #define __Pyx_PyLong_Digits(x)  (((PyLongObject*)x)->long_value.ob_digit)
  #else
  #define __Pyx_PyLong_Digits(x)  (((PyLongObject*)x)->ob_digit)
  #endif
#endif

#if __PYX_DEFAULT_STRING_ENCODING_IS_UTF8
  #define __Pyx_PyUnicode_FromStringAndSize(c_str, size) PyUnicode_DecodeUTF8(c_str, size, NULL)
#elif __PYX_DEFAULT_STRING_ENCODING_IS_ASCII
  #define __Pyx_PyUnicode_FromStringAndSize(c_str, size) PyUnicode_DecodeASCII(c_str, size, NULL)
#else
  #define __Pyx_PyUnicode_FromStringAndSize(c_str, size) PyUnicode_Decode(c_str, size, __PYX_DEFAULT_STRING_ENCODING, NULL)
#endif

/////////////// TypeConversions ///////////////

/* Type Conversion Functions */

#include <string.h>

static CYTHON_INLINE Py_ssize_t __Pyx_ssize_strlen(const char *s) {
    size_t len = strlen(s);
    if (unlikely(len > (size_t) PY_SSIZE_T_MAX)) {
        PyErr_SetString(PyExc_OverflowError, "byte string is too long");
        return -1;
    }
    return (Py_ssize_t) len;
}

static CYTHON_INLINE PyObject* __Pyx_PyUnicode_FromString(const char* c_str) {
    Py_ssize_t len = __Pyx_ssize_strlen(c_str);
    if (unlikely(len < 0)) return NULL;
    return __Pyx_PyUnicode_FromStringAndSize(c_str, len);
}

static CYTHON_INLINE PyObject* __Pyx_PyByteArray_FromString(const char* c_str) {
    Py_ssize_t len = __Pyx_ssize_strlen(c_str);
    if (unlikely(len < 0)) return NULL;
    return PyByteArray_FromStringAndSize(c_str, len);
}

// Py3.7 returns a "const char*" for unicode strings
static CYTHON_INLINE const char* __Pyx_PyObject_AsString(PyObject* o) {
    Py_ssize_t ignore;
    return __Pyx_PyObject_AsStringAndSize(o, &ignore);
}

#if __PYX_DEFAULT_STRING_ENCODING_IS_ASCII || __PYX_DEFAULT_STRING_ENCODING_IS_UTF8
static CYTHON_INLINE const char* __Pyx_PyUnicode_AsStringAndSize(PyObject* o, Py_ssize_t *length) {
    if (unlikely(__Pyx_PyUnicode_READY(o) == -1)) return NULL;
#if CYTHON_COMPILING_IN_LIMITED_API
    {
        const char* result;
        Py_ssize_t unicode_length;
        CYTHON_MAYBE_UNUSED_VAR(unicode_length); // only for __PYX_DEFAULT_STRING_ENCODING_IS_ASCII
        #if __PYX_LIMITED_VERSION_HEX < 0x030A0000
        if (unlikely(PyArg_Parse(o, "s#", &result, length) < 0)) return NULL;
        #else
        result = PyUnicode_AsUTF8AndSize(o, length);
        #endif
        #if __PYX_DEFAULT_STRING_ENCODING_IS_ASCII
        // Pre-checking "isascii" the Limited API involves making Python function calls.
        // Therefore we do a post-check instead that the lengths of the encoded and unicode
        // strings are the same (which is only true for ascii strings).
        // If this isn't true then we've already done the encoding so it's a potential
        // performance loss, but it should be better in the successful case.
        unicode_length = PyUnicode_GetLength(o);
        if (unlikely(unicode_length < 0)) return NULL;
        if (unlikely(unicode_length != *length)) {
            // raise the error
            PyUnicode_AsASCIIString(o);
            return NULL;
        }
        #endif
        return result;
    }
#else /* CYTHON_COMPILING_IN_LIMITED_API */
#if __PYX_DEFAULT_STRING_ENCODING_IS_ASCII
    if (likely(PyUnicode_IS_ASCII(o))) {
        // cached for the lifetime of the object
        *length = PyUnicode_GET_LENGTH(o);
        return PyUnicode_AsUTF8(o);
    } else {
        // raise the error
        PyUnicode_AsASCIIString(o);
        return NULL;
    }
#else
    return PyUnicode_AsUTF8AndSize(o, length);
#endif /* __PYX_DEFAULT_STRING_ENCODING_IS_ASCII */
#endif /* !CYTHON_COMPILING_IN_LIMITED_API */
}
#endif

// Py3.7 returns a "const char*" for unicode strings
static CYTHON_INLINE const char* __Pyx_PyObject_AsStringAndSize(PyObject* o, Py_ssize_t *length) {
#if __PYX_DEFAULT_STRING_ENCODING_IS_ASCII || __PYX_DEFAULT_STRING_ENCODING_IS_UTF8
    if (PyUnicode_Check(o)) {
        return __Pyx_PyUnicode_AsStringAndSize(o, length);
    } else
#endif

    if (PyByteArray_Check(o)) {
#if (CYTHON_ASSUME_SAFE_SIZE && CYTHON_ASSUME_SAFE_MACROS) || (CYTHON_COMPILING_IN_PYPY && (defined(PyByteArray_AS_STRING) && defined(PyByteArray_GET_SIZE)))
        *length = PyByteArray_GET_SIZE(o);
        return PyByteArray_AS_STRING(o);
#else
        *length = PyByteArray_Size(o);
        if (*length == -1) return NULL;
        return PyByteArray_AsString(o);
#endif

    } else
    {
        char* result;
        int r = PyBytes_AsStringAndSize(o, &result, length);
        if (unlikely(r < 0)) {
            return NULL;
        } else {
            return result;
        }
    }
}

/* Note: __Pyx_PyObject_IsTrue is written to minimize branching. */
static CYTHON_INLINE int __Pyx_PyObject_IsTrue(PyObject* x) {
   int is_true = x == Py_True;
   if (is_true | (x == Py_False) | (x == Py_None)) return is_true;
   else return PyObject_IsTrue(x);
}

static CYTHON_INLINE int __Pyx_PyObject_IsTrueAndDecref(PyObject* x) {
    int retval;
    if (unlikely(!x)) return -1;
    retval = __Pyx_PyObject_IsTrue(x);
    Py_DECREF(x);
    return retval;
}

static PyObject* __Pyx_PyNumber_LongWrongResultType(PyObject* result) {
    __Pyx_TypeName result_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(result));
    if (PyLong_Check(result)) {
        // CPython issue #17576: warn if 'result' not of exact type int.
        if (PyErr_WarnFormat(PyExc_DeprecationWarning, 1,
                "__int__ returned non-int (type " __Pyx_FMT_TYPENAME ").  "
                "The ability to return an instance of a strict subclass of int is deprecated, "
                "and may be removed in a future version of Python.",
                result_type_name)) {
            __Pyx_DECREF_TypeName(result_type_name);
            Py_DECREF(result);
            return NULL;
        }
        __Pyx_DECREF_TypeName(result_type_name);
        return result;
    }
    PyErr_Format(PyExc_TypeError,
                 "__int__ returned non-int (type " __Pyx_FMT_TYPENAME ")",
                 result_type_name);
    __Pyx_DECREF_TypeName(result_type_name);
    Py_DECREF(result);
    return NULL;
}

static CYTHON_INLINE PyObject* __Pyx_PyNumber_Long(PyObject* x) {
#if CYTHON_USE_TYPE_SLOTS
  PyNumberMethods *m;
#endif
  PyObject *res = NULL;
  if (likely(PyLong_Check(x)))
      return __Pyx_NewRef(x);
#if CYTHON_USE_TYPE_SLOTS
  m = Py_TYPE(x)->tp_as_number;
  if (likely(m && m->nb_int)) {
      res = m->nb_int(x);
  }
#else
  if (!PyBytes_CheckExact(x) && !PyUnicode_CheckExact(x)) {
      res = PyNumber_Long(x);
  }
#endif
  if (likely(res)) {
      if (unlikely(!PyLong_CheckExact(res))) {
          return __Pyx_PyNumber_LongWrongResultType(res);
      }
  }
  else if (!PyErr_Occurred()) {
      PyErr_SetString(PyExc_TypeError,
                      "an integer is required");
  }
  return res;
}

{{py: from Cython.Utility import pylong_join }}

static CYTHON_INLINE Py_ssize_t __Pyx_PyIndex_AsSsize_t(PyObject* b) {
  Py_ssize_t ival;
  PyObject *x;
  if (likely(PyLong_CheckExact(b))) {
    #if CYTHON_USE_PYLONG_INTERNALS
    // handle most common case first to avoid indirect branch and optimise branch prediction
    if (likely(__Pyx_PyLong_IsCompact(b))) {
        return __Pyx_PyLong_CompactValue(b);
    } else {
      const digit* digits = __Pyx_PyLong_Digits(b);
      const Py_ssize_t size = __Pyx_PyLong_SignedDigitCount(b);
      switch (size) {
         {{for _size in (2, 3, 4)}}
         {{for _case in (_size, -_size)}}
         case {{_case}}:
           if (8 * sizeof(Py_ssize_t) > {{_size}} * PyLong_SHIFT) {
             return {{'-' if _case < 0 else ''}}(Py_ssize_t) {{pylong_join(_size, 'digits', 'size_t')}};
           }
           break;
         {{endfor}}
         {{endfor}}
      }
    }
    #endif
    return PyLong_AsSsize_t(b);
  }
  x = PyNumber_Index(b);
  if (!x) return -1;
  ival = PyLong_AsSsize_t(x);
  Py_DECREF(x);
  return ival;
}


static CYTHON_INLINE Py_hash_t __Pyx_PyIndex_AsHash_t(PyObject* o) {
  if (sizeof(Py_hash_t) == sizeof(Py_ssize_t)) {
    return (Py_hash_t) __Pyx_PyIndex_AsSsize_t(o);
  } else {
    Py_ssize_t ival;
    PyObject *x;
    x = PyNumber_Index(o);
    if (!x) return -1;
    ival = PyLong_AsLong(x);
    Py_DECREF(x);
    return ival;
  }
}

static CYTHON_INLINE PyObject *__Pyx_Owned_Py_None(int b) {
    CYTHON_UNUSED_VAR(b);
    return __Pyx_NewRef(Py_None);
}

static CYTHON_INLINE PyObject * __Pyx_PyBool_FromLong(long b) {
  return __Pyx_NewRef(b ? Py_True: Py_False);
}


static CYTHON_INLINE PyObject * __Pyx_PyLong_FromSize_t(size_t ival) {
    return PyLong_FromSize_t(ival);
}


/////////////// pybuiltin_invalid ///////////////

static void __Pyx_PyBuiltin_Invalid(PyObject *obj, const char *type_name, const char *argname) {
    __Pyx_TypeName obj_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(obj));
    if (argname) {
        PyErr_Format(PyExc_TypeError,
            "Argument '%.200s' has incorrect type (expected %.200s, got " __Pyx_FMT_TYPENAME ")",
            argname, type_name, obj_type_name
        );
    } else {
        PyErr_Format(PyExc_TypeError,
            "Expected %.200s, got " __Pyx_FMT_TYPENAME,
            type_name, obj_type_name
        );
    }
    __Pyx_DECREF_TypeName(obj_type_name);
}


/////////////// pyfloat_simplify.proto ///////////////

static CYTHON_INLINE int __Pyx_PyFloat_FromNumber(PyObject **number_var, const char *argname, int accept_none); /*proto*/

/////////////// pyfloat_simplify ///////////////
//@requires: pybuiltin_invalid

static CYTHON_INLINE int __Pyx_PyFloat_FromNumber(PyObject **number_var, const char *argname, int accept_none) {
    // Convert any float-compatible Python number object into a Python float.
    // NOTE: This function decrefs 'number' if a conversion happens to replace the original object.
    PyObject *number = *number_var;
    if (likely((accept_none && number == Py_None) || PyFloat_CheckExact(number))) {
        return 0;
    }

    PyObject *float_object;
    if (likely(PyLong_CheckExact(number))) {
        double val;
#if CYTHON_USE_PYLONG_INTERNALS
        if (likely(__Pyx_PyLong_IsCompact(number))) {
            val = (double) __Pyx_PyLong_CompactValue(number);
        } else
#endif
        {
            val = PyLong_AsDouble(number);
            if (unlikely(val == -1.0 && PyErr_Occurred())) goto bad;
        }
        float_object = PyFloat_FromDouble(val);
    }
    else if (PyNumber_Check(number)) {
        // PyNumber_Float() also parses strings, which we must reject.
        float_object = PyNumber_Float(number);
    } else {
        __Pyx_PyBuiltin_Invalid(number, "float", argname);
        goto bad;
    }
    if (unlikely(!float_object)) goto bad;

    *number_var = float_object;
    Py_DECREF(number);
    return 0;

bad:
    *number_var = NULL;
    Py_DECREF(number);
    return -1;
}


/////////////// pyint_simplify.proto ///////////////

static CYTHON_INLINE int __Pyx_PyInt_FromNumber(PyObject **number_var, const char *argname, int accept_none); /*proto*/

/////////////// pyint_simplify ///////////////
//@requires: pybuiltin_invalid

static CYTHON_INLINE int __Pyx_PyInt_FromNumber(PyObject **number_var, const char *argname, int accept_none) {
    // Convert any int-compatible Python number object into a Python int.
    // NOTE: This function decrefs 'number' if a conversion happens to replace the original object.
    PyObject *number = *number_var;
    if (likely((accept_none && number == Py_None) || PyLong_CheckExact(number))) {
        return 0;
    }

    PyObject *int_object;
    if (likely(PyNumber_Check(number))) {
        // PyNumber_Long() also parses strings, which we must reject.
        int_object = PyNumber_Long(number);
        if (unlikely(!int_object)) goto bad;
    } else {
        __Pyx_PyBuiltin_Invalid(number, "int", argname);
        goto bad;
    }

    *number_var = int_object;
    Py_DECREF(number);
    return 0;

bad:
    *number_var = NULL;
    Py_DECREF(number);
    return -1;
}


/////////////// pynumber_float.proto ///////////////

static CYTHON_INLINE PyObject* __Pyx__PyNumber_Float(PyObject* obj); /* proto */
#define __Pyx_PyNumber_Float(x) (PyFloat_CheckExact(x) ? __Pyx_NewRef(x) : __Pyx__PyNumber_Float(x))

/////////////// pynumber_float ///////////////
//@requires: Optimize.c::pybytes_as_double
//@requires: Optimize.c::pyunicode_as_double

static CYTHON_INLINE PyObject* __Pyx__PyNumber_Float(PyObject* obj) {
    // 'obj is PyFloat' is handled in the calling macro
    double val;
    if (PyLong_CheckExact(obj)) {
#if CYTHON_USE_PYLONG_INTERNALS
        if (likely(__Pyx_PyLong_IsCompact(obj))) {
            val = (double) __Pyx_PyLong_CompactValue(obj);
            goto no_error;
        }
#endif
        val = PyLong_AsDouble(obj);
    } else if (PyUnicode_CheckExact(obj)) {
        val = __Pyx_PyUnicode_AsDouble(obj);
    } else if (PyBytes_CheckExact(obj)) {
        val = __Pyx_PyBytes_AsDouble(obj);
    } else if (PyByteArray_CheckExact(obj)) {
        val = __Pyx_PyByteArray_AsDouble(obj);
    } else {
        return PyNumber_Float(obj);
    }

    if (unlikely(val == -1 && PyErr_Occurred())) {
        return NULL;
    }
#if CYTHON_USE_PYLONG_INTERNALS
no_error:
#endif
    return PyFloat_FromDouble(val);
}

/////////////// GCCDiagnostics.proto ///////////////

// GCC diagnostic pragmas were introduced in GCC 4.6
// Used to silence conversion warnings that are ok but cannot be avoided.
#if !defined(__INTEL_COMPILER) && defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6))
#define __Pyx_HAS_GCC_DIAGNOSTIC
#endif


/////////////// ToPyCTupleUtility.proto ///////////////

static PyObject* {{funcname}}({{struct_type_decl}});

/////////////// ToPyCTupleUtility ///////////////

static PyObject* {{funcname}}({{struct_type_decl}} value) {
    PyObject* items[{{size}}] = { {{ ', '.join('0'*size) }} };
    PyObject* result = NULL;

    {{for ix, component in enumerate(components):}}
        {{py:attr = f"value.f{ix}" }}
        items[{{ix}}] = {{component.to_py_function}}({{attr}});
        if (unlikely(!items[{{ix}}])) goto bad;
    {{endfor}}

    result = PyTuple_New({{size}});
    if (unlikely(!result)) goto bad;

    for (Py_ssize_t i=0; i<{{ size }}; ++i) {
        PyObject *item = items[i];
        items[i] = NULL;
        #if !CYTHON_ASSUME_SAFE_MACROS
        // PyTuple_SetItem() always steals a reference.
        if (unlikely(PyTuple_SetItem(result, i, item) < 0)) goto bad;
        #else
        PyTuple_SET_ITEM(result, i, item);
        #endif
    }

    return result;

bad:
    Py_XDECREF(result);
    for (Py_ssize_t i={{ size-1 }}; i >= 0; --i) {
        Py_XDECREF(items[i]);
    }
    return NULL;
}


/////////////// FromPyCTupleUtility.proto ///////////////
static CYTHON_INLINE {{struct_type_decl}} {{funcname}}(PyObject *);

/////////////// FromPyCTupleUtility ///////////////

#if CYTHON_ASSUME_SAFE_MACROS && CYTHON_ASSUME_SAFE_SIZE && !CYTHON_AVOID_BORROWED_REFS
static void __Pyx_tuple_{{funcname}}(PyObject * o, {{struct_type_decl}} *result) {
    {{for ix, component in enumerate(components):}}
        {{py:attr = "result->f%s" % ix}}
        {{attr}} = {{component.from_py_function}}(PyTuple_GET_ITEM(o, {{ix}}));
        if ({{component.error_condition(attr)}}) goto bad;
    {{endfor}}
    return;
bad:
    return;
}

static void __Pyx_list_{{funcname}}(PyObject * o, {{struct_type_decl}} *result) {
    {{for ix, component in enumerate(components):}}
        {{py:attr = "result->f%s" % ix}}
        {{attr}} = {{component.from_py_function}}(PyList_GET_ITEM(o, {{ix}}));
        if ({{component.error_condition(attr)}}) goto bad;
    {{endfor}}
    return;
bad:
    return;
}
#endif

static void __Pyx_seq_{{funcname}}(PyObject * o, {{struct_type_decl}} *result) {
    if (unlikely(!PySequence_Check(o))) {
        __Pyx_TypeName o_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(o));
        PyErr_Format(PyExc_TypeError,
                     "Expected a sequence of size %zd, got " __Pyx_FMT_TYPENAME, (Py_ssize_t) {{size}}, o_type_name);
        __Pyx_DECREF_TypeName(o_type_name);
        goto bad;
    } else if (unlikely(PySequence_Length(o) != {{size}})) {
        PyErr_Format(PyExc_TypeError,
                     "Expected a sequence of size %zd, got size %zd", (Py_ssize_t) {{size}}, PySequence_Length(o));
        goto bad;
    }

    {
        PyObject *item;
    {{for ix, component in enumerate(components):}}
        {{py:attr = "result->f%s" % ix}}
        item = __Pyx_PySequence_ITEM(o, {{ix}});  if (unlikely(!item)) goto bad;
        {{attr}} = {{component.from_py_function}}(item);
        Py_DECREF(item);
        if ({{component.error_condition(attr)}}) goto bad;
    {{endfor}}
    }
    return;
bad:
    return;
}

static CYTHON_INLINE {{struct_type_decl}} {{funcname}}(PyObject * o) {
    {{struct_type_decl}} result;

    #if CYTHON_ASSUME_SAFE_MACROS && CYTHON_ASSUME_SAFE_SIZE && !CYTHON_AVOID_BORROWED_REFS
    if (likely(PyTuple_Check(o) && PyTuple_GET_SIZE(o) == {{size}})) {
        __Pyx_tuple_{{funcname}}(o, &result);
    } else if (likely(PyList_Check(o) && PyList_GET_SIZE(o) == {{size}})) {
        __Pyx_list_{{funcname}}(o, &result);
    } else
    #endif
    {
        __Pyx_seq_{{funcname}}(o, &result);
    }

    return result;
}


/////////////// UnicodeAsUCS4.proto ///////////////

static CYTHON_INLINE Py_UCS4 __Pyx_PyUnicode_AsPy_UCS4(PyObject*);

/////////////// UnicodeAsUCS4 ///////////////

static void __Pyx_PyUnicode_AsPy_UCS4_error(Py_ssize_t length) {
    if (likely(length >= 0)) {
        // "length == -1" indicates an error already.
        PyErr_Format(PyExc_ValueError,
                     "only single character unicode strings can be converted to Py_UCS4, "
                     "got length %" CYTHON_FORMAT_SSIZE_T "d", length);
    }
}

static CYTHON_INLINE Py_UCS4 __Pyx_PyUnicode_AsPy_UCS4(PyObject* x) {
    Py_ssize_t length = __Pyx_PyUnicode_GET_LENGTH(x);
    if (unlikely(length != 1)) {
        __Pyx_PyUnicode_AsPy_UCS4_error(length);
        return (Py_UCS4)-1;
    }
    return __Pyx_PyUnicode_READ_CHAR(x, 0);
}


/////////////// ObjectAsUCS4.proto ///////////////
//@requires: UnicodeAsUCS4

static Py_UCS4 __Pyx__PyObject_AsPy_UCS4(PyObject*);

static CYTHON_INLINE Py_UCS4 __Pyx_PyObject_AsPy_UCS4(PyObject *x) {
    return (likely(PyUnicode_Check(x)) ? __Pyx_PyUnicode_AsPy_UCS4(x) : __Pyx__PyObject_AsPy_UCS4(x));
}


/////////////// ObjectAsUCS4 ///////////////

static void __Pyx__PyObject_AsPy_UCS4_raise_error(long ival) {
   if (ival < 0) {
       if (!PyErr_Occurred())
           PyErr_SetString(PyExc_OverflowError,
                           "cannot convert negative value to Py_UCS4");
   } else {
       PyErr_SetString(PyExc_OverflowError,
                       "value too large to convert to Py_UCS4");
   }
}

static Py_UCS4 __Pyx__PyObject_AsPy_UCS4(PyObject* x) {
   long ival;
   ival = __Pyx_PyLong_As_long(x);
   if (unlikely(!__Pyx_is_valid_index(ival, 1114111 + 1))) {
       __Pyx__PyObject_AsPy_UCS4_raise_error(ival);
       return (Py_UCS4)-1;
   }
   return (Py_UCS4)ival;
}


/////////////// ObjectAsPyUnicode.proto ///////////////

static CYTHON_INLINE Py_UNICODE __Pyx_PyObject_AsPy_UNICODE(PyObject*);

/////////////// ObjectAsPyUnicode ///////////////

static CYTHON_INLINE Py_UNICODE __Pyx_PyObject_AsPy_UNICODE(PyObject* x) {
    long ival;
  #if defined(Py_UNICODE_SIZE) && Py_UNICODE_SIZE == 2
    const long maxval = 65535;
  #else
    const long maxval = 1114111;
  #endif
    if (PyUnicode_Check(x)) {
        Py_ssize_t length = __Pyx_PyUnicode_GET_LENGTH(x);
        if (unlikely(length != 1)) {
            // -1 indicates an error.
            if (length >= 0) {
                PyErr_Format(PyExc_ValueError,
                             "only single character unicode strings can be converted to Py_UNICODE, "
                             "got length %" CYTHON_FORMAT_SSIZE_T "d", length);
            }
            return (Py_UNICODE)-1;
        }
        ival = PyUnicode_READ_CHAR(x, 0);
    } else {
        ival = __Pyx_PyLong_As_long(x);
    }
    if (unlikely(!__Pyx_is_valid_index(ival, maxval + 1))) {
        if (ival < 0) {
            if (!PyErr_Occurred())
                PyErr_SetString(PyExc_OverflowError,
                                "cannot convert negative value to Py_UNICODE");
            return (Py_UNICODE)-1;
        } else {
            PyErr_SetString(PyExc_OverflowError,
                            "value too large to convert to Py_UNICODE");
        }
        return (Py_UNICODE)-1;
    }
    return (Py_UNICODE)ival;
}


/////////////// CIntToPy.proto ///////////////

static CYTHON_INLINE PyObject* {{TO_PY_FUNCTION}}({{TYPE}} value);

/////////////// CIntToPy ///////////////
//@requires: GCCDiagnostics
//@requires: ObjectHandling.c::PyObjectVectorCallKwBuilder

static CYTHON_INLINE PyObject* {{TO_PY_FUNCTION}}({{TYPE}} value) {
#ifdef __Pyx_HAS_GCC_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif
    const {{TYPE}} neg_one = ({{TYPE}}) -1, const_zero = ({{TYPE}}) 0;
#ifdef __Pyx_HAS_GCC_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
    const int is_unsigned = neg_one > const_zero;
    if (is_unsigned) {
        if (sizeof({{TYPE}}) < sizeof(long)) {
            return PyLong_FromLong((long) value);
        } else if (sizeof({{TYPE}}) <= sizeof(unsigned long)) {
            return PyLong_FromUnsignedLong((unsigned long) value);
#if !CYTHON_COMPILING_IN_PYPY
        // PyLong_FromUnsignedLongLong() does not necessarily accept ULL arguments in PyPy.
        // See https://github.com/cython/cython/issues/6890
        } else if (sizeof({{TYPE}}) <= sizeof(unsigned PY_LONG_LONG)) {
            return PyLong_FromUnsignedLongLong((unsigned PY_LONG_LONG) value);
#endif
        }
    } else {
        if (sizeof({{TYPE}}) <= sizeof(long)) {
            return PyLong_FromLong((long) value);
        } else if (sizeof({{TYPE}}) <= sizeof(PY_LONG_LONG)) {
            return PyLong_FromLongLong((PY_LONG_LONG) value);
        }
    }
    {
        unsigned char *bytes = (unsigned char *)&value;
#if !CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX >= 0x030d00A4
        if (is_unsigned) {
            return PyLong_FromUnsignedNativeBytes(bytes, sizeof(value), -1);
        } else {
            return PyLong_FromNativeBytes(bytes, sizeof(value), -1);
        }
#elif !CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX < 0x030d0000
        int one = 1; int little = (int)*(unsigned char *)&one;
        return _PyLong_FromByteArray(bytes, sizeof({{TYPE}}),
                                     little, !is_unsigned);
#else
        // call int.from_bytes()
        int one = 1; int little = (int)*(unsigned char *)&one;
        PyObject *from_bytes, *result = NULL, *kwds = NULL;
        PyObject *py_bytes = NULL, *order_str = NULL;
        from_bytes = PyObject_GetAttrString((PyObject*)&PyLong_Type, "from_bytes");
        if (!from_bytes) return NULL;
        py_bytes = PyBytes_FromStringAndSize((char*)bytes, sizeof({{TYPE}}));
        if (!py_bytes) goto limited_bad;
        // I'm deliberately not using PYIDENT here because this code path is very unlikely
        // to ever run so it seems a pessimization mostly.
        order_str = PyUnicode_FromString(little ? "little" : "big");
        if (!order_str) goto limited_bad;
        {
            PyObject *args[3+(CYTHON_VECTORCALL ? 1 : 0)] = { NULL, py_bytes, order_str };
            if (!is_unsigned) {
                kwds = __Pyx_MakeVectorcallBuilderKwds(1);
                if (!kwds) goto limited_bad;
                if (__Pyx_VectorcallBuilder_AddArgStr("signed", __Pyx_NewRef(Py_True), kwds, args+3, 0) < 0) goto limited_bad;
            }
            result = __Pyx_Object_Vectorcall_CallFromBuilder(from_bytes, args+1, 2 | __Pyx_PY_VECTORCALL_ARGUMENTS_OFFSET, kwds);
        }

        limited_bad:
        Py_XDECREF(kwds);
        Py_XDECREF(order_str);
        Py_XDECREF(py_bytes);
        Py_XDECREF(from_bytes);
        return result;
#endif
    }
}


/////////////// COrdinalToPyUnicode.proto ///////////////

static CYTHON_INLINE int __Pyx_CheckUnicodeValue(int value);
static CYTHON_INLINE PyObject* __Pyx_PyUnicode_FromOrdinal_Padded(int value, Py_ssize_t width, char padding_char);

/////////////// COrdinalToPyUnicode ///////////////
//@requires: StringTools.c::BuildPyUnicode

static CYTHON_INLINE int __Pyx_CheckUnicodeValue(int value) {
    return value <= 1114111;
}

static PyObject* __Pyx_PyUnicode_FromOrdinal_Padded(int value, Py_ssize_t ulength, char padding_char) {
    Py_ssize_t padding_length = ulength - 1;
    if (likely((padding_length <= 250) && (value < 0xD800 || value > 0xDFFF))) {
        // Encode to UTF-8 / Latin1 buffer, then decode.
        char chars[256];

        if (value <= 255) {
            // Simple Latin1 result, fast to decode.
            memset(chars, padding_char, (size_t) padding_length);
            chars[ulength-1] = (char) value;
            return PyUnicode_DecodeLatin1(chars, ulength, NULL);
        }

        char *cpos = chars + sizeof(chars);
        if (value < 0x800) {
            *--cpos = (char) (0x80 | (value & 0x3f));
            value >>= 6;
            *--cpos = (char) (0xc0 | (value & 0x1f));
        } else if (value < 0x10000) {
            *--cpos = (char) (0x80 | (value & 0x3f));
            value >>= 6;
            *--cpos = (char) (0x80 | (value & 0x3f));
            value >>= 6;
            *--cpos = (char) (0xe0 | (value & 0x0f));
        } else {
            *--cpos = (char) (0x80 | (value & 0x3f));
            value >>= 6;
            *--cpos = (char) (0x80 | (value & 0x3f));
            value >>= 6;
            *--cpos = (char) (0x80 | (value & 0x3f));
            value >>= 6;
            *--cpos = (char) (0xf0 | (value & 0x07));
        }
        cpos -= padding_length;
        memset(cpos, padding_char, (size_t) padding_length);
        return PyUnicode_DecodeUTF8(cpos, chars + sizeof(chars) - cpos, NULL);
    }

    if (value <= 127 && CYTHON_USE_UNICODE_INTERNALS) {
        const char chars[1] = {(char) value};
        return __Pyx_PyUnicode_BuildFromAscii(ulength, chars, 1, 0, padding_char);
    }

    {
        PyObject *uchar, *padding_uchar, *padding, *result;

        padding_uchar = PyUnicode_FromOrdinal(padding_char);
        if (unlikely(!padding_uchar)) return NULL;
        padding = PySequence_Repeat(padding_uchar, padding_length);
        Py_DECREF(padding_uchar);
        if (unlikely(!padding)) return NULL;

        uchar = PyUnicode_FromOrdinal(value);
        if (unlikely(!uchar)) {
            Py_DECREF(padding);
            return NULL;
        }

        result = PyUnicode_Concat(padding, uchar);
        Py_DECREF(padding);
        Py_DECREF(uchar);
        return result;
    }
}


/////////////// CIntToDigits ///////////////

static const char DIGIT_PAIRS_10[2*10*10+1] = {
    "00010203040506070809"
    "10111213141516171819"
    "20212223242526272829"
    "30313233343536373839"
    "40414243444546474849"
    "50515253545556575859"
    "60616263646566676869"
    "70717273747576777879"
    "80818283848586878889"
    "90919293949596979899"
};

static const char DIGIT_PAIRS_8[2*8*8+1] = {
    "0001020304050607"
    "1011121314151617"
    "2021222324252627"
    "3031323334353637"
    "4041424344454647"
    "5051525354555657"
    "6061626364656667"
    "7071727374757677"
};

static const char DIGITS_HEX[2*16+1] = {
    "0123456789abcdef"
    "0123456789ABCDEF"
};


/////////////// CIntToPyUnicode.proto ///////////////

#define {{TO_PY_FUNCTION}}(value, width, padding_char, format_char) ( \
    ((format_char) == ('c')) ? \
        __Pyx_uchar_{{TO_PY_FUNCTION}}(value, width, padding_char) : \
        __Pyx__{{TO_PY_FUNCTION}}(value, width, padding_char, format_char) \
    )
static CYTHON_INLINE PyObject* __Pyx_uchar_{{TO_PY_FUNCTION}}({{TYPE}} value, Py_ssize_t width, char padding_char);
static CYTHON_INLINE PyObject* __Pyx__{{TO_PY_FUNCTION}}({{TYPE}} value, Py_ssize_t width, char padding_char, char format_char);

/////////////// CIntToPyUnicode ///////////////
//@requires: StringTools.c::IncludeStringH
//@requires: StringTools.c::BuildPyUnicode
//@requires: ModuleSetupCode.c::IncludeStdlibH
//@requires: COrdinalToPyUnicode
//@requires: CIntToDigits
//@requires: GCCDiagnostics

// NOTE: inlining because most arguments are constant, which collapses lots of code below

static CYTHON_INLINE PyObject* __Pyx_uchar_{{TO_PY_FUNCTION}}({{TYPE}} value, Py_ssize_t width, char padding_char) {
#ifdef __Pyx_HAS_GCC_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif
    const {{TYPE}} neg_one = ({{TYPE}}) -1, const_zero = ({{TYPE}}) 0;
#ifdef __Pyx_HAS_GCC_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
    const int is_unsigned = neg_one > const_zero;

    // This check is just an awful variation on "(0 <= value <= 1114111)",
    // but without C compiler complaints about compile time constant conditions depending on the signed/unsigned TYPE.
    if (unlikely(!(is_unsigned || value == 0 || value > 0) ||
                    !(sizeof(value) <= 2 || value & ~ ({{TYPE}}) 0x01fffff || __Pyx_CheckUnicodeValue((int) value)))) {
        // PyUnicode_FromOrdinal() and chr() raise ValueError, f-strings raise OverflowError. :-/
        PyErr_SetString(PyExc_OverflowError, "%c arg not in range(0x110000)");
        return NULL;
    }
    if (width <= 1) {
        return PyUnicode_FromOrdinal((int) value);
    }
    return __Pyx_PyUnicode_FromOrdinal_Padded((int) value, width, padding_char);
}

static CYTHON_INLINE PyObject* __Pyx__{{TO_PY_FUNCTION}}({{TYPE}} value, Py_ssize_t width, char padding_char, char format_char) {
    // simple and conservative C string allocation on the stack: each byte gives at most 3 digits, plus sign
    char digits[sizeof({{TYPE}})*3+2];
    // 'dpos' points to end of digits array + 1 initially to allow for pre-decrement looping
    char *dpos, *end = digits + sizeof({{TYPE}})*3+2;
    const char *hex_digits = DIGITS_HEX;
    Py_ssize_t length, ulength;
    int prepend_sign, last_one_off;
    {{TYPE}} remaining;
#ifdef __Pyx_HAS_GCC_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif
    const {{TYPE}} neg_one = ({{TYPE}}) -1, const_zero = ({{TYPE}}) 0;
#ifdef __Pyx_HAS_GCC_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
    const int is_unsigned = neg_one > const_zero;

    if (format_char == 'X') {
        hex_digits += 16;
        format_char = 'x';
    }

    // surprise: even trivial sprintf() calls don't get optimised in gcc (4.8)
    remaining = value; /* not using abs(value) to avoid overflow problems */
    last_one_off = 0;
    dpos = end;
    do {
        int digit_pos;
        switch (format_char) {
        case 'o':
            digit_pos = abs((int)(remaining % (8*8)));
            remaining = ({{TYPE}}) (remaining / (8*8));
            dpos -= 2;
            memcpy(dpos, DIGIT_PAIRS_8 + digit_pos * 2, 2); /* copy 2 digits at a time, unaligned */
            last_one_off = (digit_pos < 8);
            break;
        case 'd':
            digit_pos = abs((int)(remaining % (10*10)));
            remaining = ({{TYPE}}) (remaining / (10*10));
            dpos -= 2;
            memcpy(dpos, DIGIT_PAIRS_10 + digit_pos * 2, 2); /* copy 2 digits at a time, unaligned */
            last_one_off = (digit_pos < 10);
            break;
        case 'x':
            *(--dpos) = hex_digits[abs((int)(remaining % 16))];
            remaining = ({{TYPE}}) (remaining / 16);
            break;
        default:
            assert(0);
            break;
        }
    } while (unlikely(remaining != 0));

    // Correct dpos by 1 if we read an excess digit.
    assert(!last_one_off || *dpos == '0');
    dpos += last_one_off;

    length = end - dpos;
    ulength = length;
    prepend_sign = 0;
    if (!is_unsigned && value <= neg_one) {
        if (padding_char == ' ' || width <= length + 1) {
            *(--dpos) = '-';
            ++length;
        } else {
            prepend_sign = 1;
        }
        ++ulength;
    }
    if (width > ulength) {
        ulength = width;
    }
    // single character unicode strings are cached in CPython => use PyUnicode_FromOrdinal() for them
    if (ulength == 1) {
        return PyUnicode_FromOrdinal(*dpos);
    }
    return __Pyx_PyUnicode_BuildFromAscii(ulength, dpos, (int) length, prepend_sign, padding_char);
}


/////////////// CBIntToPyUnicode.proto ///////////////

#define {{TO_PY_FUNCTION}}(value)  \
    ((value) ? __Pyx_NewRef({{TRUE_CONST}}) : __Pyx_NewRef({{FALSE_CONST}}))


/////////////// CIntFromPyVerify ///////////////

// see CIntFromPy
#define __PYX_VERIFY_RETURN_INT(target_type, func_type, func_value)       \
    __PYX__VERIFY_RETURN_INT(target_type, func_type, func_value, 0)

#define __PYX_VERIFY_RETURN_INT_EXC(target_type, func_type, func_value)   \
    __PYX__VERIFY_RETURN_INT(target_type, func_type, func_value, 1)

#define __PYX__VERIFY_RETURN_INT(target_type, func_type, func_value, exc) \
    {                                                                     \
        func_type value = func_value;                                     \
        if (sizeof(target_type) < sizeof(func_type)) {                    \
            if (unlikely(value != (func_type) (target_type) value)) {     \
                func_type zero = 0;                                       \
                if (exc && unlikely(value == (func_type)-1 && PyErr_Occurred()))  \
                    return (target_type) -1;                              \
                if (is_unsigned && unlikely(value < zero))                \
                    goto raise_neg_overflow;                              \
                else                                                      \
                    goto raise_overflow;                                  \
            }                                                             \
        }                                                                 \
        return (target_type) value;                                       \
    }


/////////////// CIntFromPy.proto ///////////////

static CYTHON_INLINE {{TYPE}} {{FROM_PY_FUNCTION}}(PyObject *);

/////////////// CIntFromPy ///////////////
//@requires: CIntFromPyVerify
//@requires: GCCDiagnostics

{{py: from Cython.Utility import pylong_join }}

static CYTHON_INLINE {{TYPE}} {{FROM_PY_FUNCTION}}(PyObject *x) {
#ifdef __Pyx_HAS_GCC_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif
    const {{TYPE}} neg_one = ({{TYPE}}) -1, const_zero = ({{TYPE}}) 0;
#ifdef __Pyx_HAS_GCC_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
    const int is_unsigned = neg_one > const_zero;

    if (unlikely(!PyLong_Check(x))) {
        {{TYPE}} val;
        PyObject *tmp = __Pyx_PyNumber_Long(x);
        if (!tmp) return ({{TYPE}}) -1;
        val = {{FROM_PY_FUNCTION}}(tmp);
        Py_DECREF(tmp);
        return val;
    }

    if (is_unsigned) {
#if CYTHON_USE_PYLONG_INTERNALS
        if (unlikely(__Pyx_PyLong_IsNeg(x))) {
            goto raise_neg_overflow;
        //} else if (__Pyx_PyLong_IsZero(x)) {
        //    return ({{TYPE}}) 0;
        } else if (__Pyx_PyLong_IsCompact(x)) {
            __PYX_VERIFY_RETURN_INT({{TYPE}}, __Pyx_compact_upylong, __Pyx_PyLong_CompactValueUnsigned(x))
        } else {
            const digit* digits = __Pyx_PyLong_Digits(x);
            assert(__Pyx_PyLong_DigitCount(x) > 1);
            switch (__Pyx_PyLong_DigitCount(x)) {
                {{for _size in (2, 3, 4)}}
                case {{_size}}:
                    if ((8 * sizeof({{TYPE}}) > {{_size-1}} * PyLong_SHIFT)) {
                        if ((8 * sizeof(unsigned long) > {{_size}} * PyLong_SHIFT)) {
                            __PYX_VERIFY_RETURN_INT({{TYPE}}, unsigned long, {{pylong_join(_size, 'digits')}})
                        } else if ((8 * sizeof({{TYPE}}) >= {{_size}} * PyLong_SHIFT)) {
                            return ({{TYPE}}) {{pylong_join(_size, 'digits', TYPE)}};
                        }
                    }
                    break;
                {{endfor}}
            }
        }
#endif
#if CYTHON_COMPILING_IN_CPYTHON && PY_VERSION_HEX < 0x030C00A7
        if (unlikely(Py_SIZE(x) < 0)) {
            goto raise_neg_overflow;
        }
#else
        {
            // misuse Py_False as a quick way to compare to a '0' int object in PyPy
            int result = PyObject_RichCompareBool(x, Py_False, Py_LT);
            if (unlikely(result < 0))
                return ({{TYPE}}) -1;
            if (unlikely(result == 1))
                goto raise_neg_overflow;
        }
#endif
        if ((sizeof({{TYPE}}) <= sizeof(unsigned long))) {
            __PYX_VERIFY_RETURN_INT_EXC({{TYPE}}, unsigned long, PyLong_AsUnsignedLong(x))
        } else if ((sizeof({{TYPE}}) <= sizeof(unsigned PY_LONG_LONG))) {
            __PYX_VERIFY_RETURN_INT_EXC({{TYPE}}, unsigned PY_LONG_LONG, PyLong_AsUnsignedLongLong(x))
        }

    } else {
        // signed
#if CYTHON_USE_PYLONG_INTERNALS
        if (__Pyx_PyLong_IsCompact(x)) {
            __PYX_VERIFY_RETURN_INT({{TYPE}}, __Pyx_compact_pylong, __Pyx_PyLong_CompactValue(x))
        } else {
            const digit* digits = __Pyx_PyLong_Digits(x);
            assert(__Pyx_PyLong_DigitCount(x) > 1);
            switch (__Pyx_PyLong_SignedDigitCount(x)) {
                {{for _size in (2, 3, 4)}}
                {{for _case in (-_size, _size)}}
                case {{_case}}:
                    if ((8 * sizeof({{TYPE}}){{' - 1' if _case < 0 else ''}} > {{_size-1}} * PyLong_SHIFT)) {
                        if ((8 * sizeof(unsigned long) > {{_size}} * PyLong_SHIFT)) {
                            __PYX_VERIFY_RETURN_INT({{TYPE}}, {{'long' if _case < 0 else 'unsigned long'}}, {{'-(long) ' if _case < 0 else ''}}{{pylong_join(_size, 'digits')}})
                        } else if ((8 * sizeof({{TYPE}}) - 1 > {{_size}} * PyLong_SHIFT)) {
                            return ({{TYPE}}) ({{'((%s)-1)*' % TYPE if _case < 0 else ''}}{{pylong_join(_size, 'digits', TYPE)}});
                        }
                    }
                    break;
                {{endfor}}
                {{endfor}}
            }
        }
#endif
        if ((sizeof({{TYPE}}) <= sizeof(long))) {
            __PYX_VERIFY_RETURN_INT_EXC({{TYPE}}, long, PyLong_AsLong(x))
        } else if ((sizeof({{TYPE}}) <= sizeof(PY_LONG_LONG))) {
            __PYX_VERIFY_RETURN_INT_EXC({{TYPE}}, PY_LONG_LONG, PyLong_AsLongLong(x))
        }
    }

    // large integer type and no access to PyLong internals => allow for a more expensive conversion
    {
        {{TYPE}} val;
        int ret = -1;
#if PY_VERSION_HEX >= 0x030d00A6 && !CYTHON_COMPILING_IN_LIMITED_API
        Py_ssize_t bytes_copied = PyLong_AsNativeBytes(
            x, &val, sizeof(val), Py_ASNATIVEBYTES_NATIVE_ENDIAN | (is_unsigned ? Py_ASNATIVEBYTES_UNSIGNED_BUFFER | Py_ASNATIVEBYTES_REJECT_NEGATIVE : 0));
        if (unlikely(bytes_copied == -1)) {
            // failed
        } else if (unlikely(bytes_copied > (Py_ssize_t) sizeof(val))) {
            goto raise_overflow;
        } else {
            ret = 0;
        }
#elif PY_VERSION_HEX < 0x030d0000 && !(CYTHON_COMPILING_IN_PYPY || CYTHON_COMPILING_IN_LIMITED_API) || defined(_PyLong_AsByteArray)
        int one = 1; int is_little = (int)*(unsigned char *)&one;
        unsigned char *bytes = (unsigned char *)&val;
        ret = _PyLong_AsByteArray((PyLongObject *)x,
                                    bytes, sizeof(val),
                                    is_little, !is_unsigned);
#else
{{if IS_ENUM}}
        // The fallback implementation uses math operations like shifting, which do not work well with enums.
        PyErr_SetString(PyExc_RuntimeError,
                        "_PyLong_AsByteArray() or PyLong_AsNativeBytes() not available, cannot convert large enums");
        val = ({{TYPE}}) -1;
{{else}}
// Inefficient copy of bit chunks through the C-API.  Probably still better than a "cannot do this" exception.
// This is substantially faster in CPython (>30%) than calling "int.to_bytes()" through the C-API.
        PyObject *v;
        PyObject *stepval = NULL, *mask = NULL, *shift = NULL;
        int bits, remaining_bits, is_negative = 0;
        int chunk_size = (sizeof(long) < 8) ? 30 : 62;

        // Use exact PyLong to prevent user defined &&/<</etc. implementations (and make Py_SIZE() work below).
        if (likely(PyLong_CheckExact(x))) {
            v = __Pyx_NewRef(x);
        } else {
            v = PyNumber_Long(x);
            if (unlikely(!v)) return ({{TYPE}}) -1;
            assert(PyLong_CheckExact(v));
        }

        // Misuse Py_False as a quick way to compare to a '0' int object.
        {
            int result = PyObject_RichCompareBool(v, Py_False, Py_LT);
            if (unlikely(result < 0)) {
                Py_DECREF(v);
                return ({{TYPE}}) -1;
            }
            is_negative = result == 1;
        }

        if (is_unsigned && unlikely(is_negative)) {
            Py_DECREF(v);
            goto raise_neg_overflow;
        } else if (is_negative) {
            // bit-invert to make sure we can safely convert it
            stepval = PyNumber_Invert(v);
            Py_DECREF(v);
            if (unlikely(!stepval))
                return ({{TYPE}}) -1;
        } else {
            stepval = v;
        }
        v = NULL;

        // Unpack full chunks of bits.
        val = ({{TYPE}}) 0;
        mask = PyLong_FromLong((1L << chunk_size) - 1); if (unlikely(!mask)) goto done;
        shift = PyLong_FromLong(chunk_size); if (unlikely(!shift)) goto done;
        for (bits = 0; bits < (int) sizeof({{TYPE}}) * 8 - chunk_size; bits += chunk_size) {
            PyObject *tmp, *digit;
            long idigit;

            digit = PyNumber_And(stepval, mask);
            if (unlikely(!digit)) goto done;

            idigit = PyLong_AsLong(digit);
            Py_DECREF(digit);
            if (unlikely(idigit < 0)) goto done;
            val |= (({{TYPE}}) idigit) << bits;

            tmp = PyNumber_Rshift(stepval, shift);
            if (unlikely(!tmp)) goto done;
            Py_DECREF(stepval); stepval = tmp;
        }

        Py_DECREF(shift); shift = NULL;
        Py_DECREF(mask); mask = NULL;

        // Add the last bits and detect overflow.
        {
            long idigit = PyLong_AsLong(stepval);
            if (unlikely(idigit < 0)) goto done;
            remaining_bits = ((int) sizeof({{TYPE}}) * 8) - bits - (is_unsigned ? 0 : 1);
            if (unlikely(idigit >= (1L << remaining_bits)))
                goto raise_overflow;
            val |= (({{TYPE}}) idigit) << bits;
        }

        // Handle sign and overflow into sign bit.
        if (!is_unsigned) {
            // gcc warns about unsigned (val < 0) => test sign bit instead
            if (unlikely(val & ((({{TYPE}}) 1) << (sizeof({{TYPE}}) * 8 - 1))))
                goto raise_overflow;
            // undo the PyNumber_Invert() above
            if (is_negative)
                val = ~val;
        }
        ret = 0;
    done:
        Py_XDECREF(shift);
        Py_XDECREF(mask);
        Py_XDECREF(stepval);
{{endif}}
#endif
        if (unlikely(ret))
            return ({{TYPE}}) -1;
        return val;
    }

raise_overflow:
    PyErr_SetString(PyExc_OverflowError,
        "value too large to convert to {{TYPE}}");
    return ({{TYPE}}) -1;

raise_neg_overflow:
    PyErr_SetString(PyExc_OverflowError,
        "can't convert negative value to {{TYPE}}");
    return ({{TYPE}}) -1;
}
