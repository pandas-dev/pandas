/*
 * Special implementations of built-in functions and methods.
 *
 * Optional optimisations for builtins are in Optimize.c.
 *
 * General object operations and protocols are in ObjectHandling.c.
 */

//////////////////// Globals.proto ////////////////////

static PyObject* __Pyx_Globals(void); /*proto*/

//////////////////// Globals ////////////////////
//@requires: ObjectHandling.c::GetAttr

// This is a stub implementation until we have something more complete.
// Currently, we only handle the most common case of a read-only dict
// of Python names.  Supporting cdef names in the module and write
// access requires a rewrite as a dedicated class.

static PyObject* __Pyx_Globals(void) {
    return __Pyx_NewRef(NAMED_CGLOBAL(moddict_cname));
}

//////////////////// PyExecGlobals.proto ////////////////////

static PyObject* __Pyx_PyExecGlobals(PyObject*);

//////////////////// PyExecGlobals ////////////////////
//@requires: PyExec

static PyObject* __Pyx_PyExecGlobals(PyObject* code) {
    return __Pyx_PyExec2(code, NAMED_CGLOBAL(moddict_cname));
}

//////////////////// PyExec.proto ////////////////////

static PyObject* __Pyx_PyExec3(PyObject*, PyObject*, PyObject*);
static CYTHON_INLINE PyObject* __Pyx_PyExec2(PyObject*, PyObject*);

//////////////////// PyExec ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyExec2(PyObject* o, PyObject* globals) {
    return __Pyx_PyExec3(o, globals, NULL);
}

static PyObject* __Pyx_PyExec3(PyObject* o, PyObject* globals, PyObject* locals) {
    PyObject* result;
#if !CYTHON_COMPILING_IN_LIMITED_API
    PyObject* s = 0;
    char *code = 0;
#endif

    if (!globals || globals == Py_None) {
        globals = NAMED_CGLOBAL(moddict_cname);
    }
#if !CYTHON_COMPILING_IN_LIMITED_API
    // In Limited API we just use exec builtin which already has this
    else if (unlikely(!PyDict_Check(globals))) {
        __Pyx_TypeName globals_type_name =
            __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(globals));
        PyErr_Format(PyExc_TypeError,
                     "exec() arg 2 must be a dict, not " __Pyx_FMT_TYPENAME,
                     globals_type_name);
        __Pyx_DECREF_TypeName(globals_type_name);
        goto bad;
    }
#endif
    if (!locals || locals == Py_None) {
        locals = globals;
    }

#if !CYTHON_COMPILING_IN_LIMITED_API
    if (__Pyx_PyDict_GetItemStr(globals, PYIDENT("__builtins__")) == NULL) {
        if (unlikely(PyDict_SetItem(globals, PYIDENT("__builtins__"), PyEval_GetBuiltins()) < 0))
            goto bad;
    }

    if (PyCode_Check(o)) {
        if (unlikely(__Pyx_PyCode_HasFreeVars((PyCodeObject *)o))) {
            PyErr_SetString(PyExc_TypeError,
                "code object passed to exec() may not contain free variables");
            goto bad;
        }
        #if CYTHON_COMPILING_IN_PYPY && PYPY_VERSION_NUM < 0x07030400
        result = PyEval_EvalCode((PyCodeObject *)o, globals, locals);
        #else
        result = PyEval_EvalCode(o, globals, locals);
        #endif
    } else {
        PyCompilerFlags cf;
        cf.cf_flags = 0;
        cf.cf_feature_version = PY_MINOR_VERSION;
        if (PyUnicode_Check(o)) {
            cf.cf_flags = PyCF_SOURCE_IS_UTF8;
            s = PyUnicode_AsUTF8String(o);
            if (unlikely(!s)) goto bad;
            o = s;
        } else if (unlikely(!PyBytes_Check(o))) {
            __Pyx_TypeName o_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(o));
            PyErr_Format(PyExc_TypeError,
                "exec: arg 1 must be string, bytes or code object, got " __Pyx_FMT_TYPENAME,
                o_type_name);
            __Pyx_DECREF_TypeName(o_type_name);
            goto bad;
        }
        code = PyBytes_AS_STRING(o);
        if (PyEval_MergeCompilerFlags(&cf)) {
            result = PyRun_StringFlags(code, Py_file_input, globals, locals, &cf);
        } else {
            result = PyRun_String(code, Py_file_input, globals, locals);
        }
        Py_XDECREF(s);
    }

    return result;
bad:
    Py_XDECREF(s);
    return 0;
#else // CYTHON_COMPILING_IN_LIMITED_API
    {
        // For the limited API we just defer to the actual builtin
        // (after setting up globals and locals) - there's too much we can't do otherwise
        PyObject *builtins, *exec, *exec_str;
        builtins = PyEval_GetBuiltins();
        if (!builtins) return NULL;
        exec_str = PyUnicode_FromStringAndSize("exec", 4);
        if (!exec_str) return NULL;
        exec = PyObject_GetItem(builtins, exec_str);
        Py_DECREF(exec_str);
        if (!exec) return NULL;
        result = PyObject_CallFunctionObjArgs(exec, o, globals, locals, NULL);
        Py_DECREF(exec);
        return result;
    }
#endif
}

//////////////////// GetAttr3.proto ////////////////////

static CYTHON_INLINE PyObject *__Pyx_GetAttr3(PyObject *, PyObject *, PyObject *); /*proto*/

//////////////////// GetAttr3 ////////////////////
//@requires: ObjectHandling.c::PyObjectGetAttrStr
//@requires: Exceptions.c::PyThreadStateGet
//@requires: Exceptions.c::PyErrFetchRestore
//@requires: Exceptions.c::PyErrExceptionMatches

#if __PYX_LIMITED_VERSION_HEX < 0x030d0000
static PyObject *__Pyx_GetAttr3Default(PyObject *d) {
    __Pyx_PyThreadState_declare
    __Pyx_PyThreadState_assign
    if (unlikely(!__Pyx_PyErr_ExceptionMatches(PyExc_AttributeError)))
        return NULL;
    __Pyx_PyErr_Clear();
    Py_INCREF(d);
    return d;
}
#endif

static CYTHON_INLINE PyObject *__Pyx_GetAttr3(PyObject *o, PyObject *n, PyObject *d) {
    PyObject *r;
#if __PYX_LIMITED_VERSION_HEX >= 0x030d0000
    int res = PyObject_GetOptionalAttr(o, n, &r);
    // On failure (res == -1), r is set to NULL.
    return (res != 0) ? r : __Pyx_NewRef(d);
#else
  #if CYTHON_USE_TYPE_SLOTS
    if (likely(PyUnicode_Check(n))) {
        r = __Pyx_PyObject_GetAttrStrNoError(o, n);
        if (unlikely(!r) && likely(!PyErr_Occurred())) {
            r = __Pyx_NewRef(d);
        }
        return r;
    }
  #endif
    r = PyObject_GetAttr(o, n);
    return (likely(r)) ? r : __Pyx_GetAttr3Default(d);
#endif
}

//////////////////// HasAttr.proto ////////////////////

#if __PYX_LIMITED_VERSION_HEX >= 0x030d0000
#define __Pyx_HasAttr(o, n)  PyObject_HasAttrWithError(o, n)
#else
static CYTHON_INLINE int __Pyx_HasAttr(PyObject *, PyObject *); /*proto*/
#endif

//////////////////// HasAttr ////////////////////
//@requires: ObjectHandling.c::PyObjectGetAttrStrNoError

#if __PYX_LIMITED_VERSION_HEX < 0x030d0000
static CYTHON_INLINE int __Pyx_HasAttr(PyObject *o, PyObject *n) {
    PyObject *r;
    if (unlikely(!PyUnicode_Check(n))) {
        PyErr_SetString(PyExc_TypeError,
                        "hasattr(): attribute name must be string");
        return -1;
    }
    r = __Pyx_PyObject_GetAttrStrNoError(o, n);
    if (!r) {
        return (unlikely(PyErr_Occurred())) ? -1 : 0;
    } else {
        Py_DECREF(r);
        return 1;
    }
}
#endif

//////////////////// Intern.proto ////////////////////

static PyObject* __Pyx_Intern(PyObject* s); /* proto */

//////////////////// Intern ////////////////////
//@requires: ObjectHandling.c::RaiseUnexpectedTypeError

static PyObject* __Pyx_Intern(PyObject* s) {
    if (unlikely(!PyUnicode_CheckExact(s))) {
        __Pyx_RaiseUnexpectedTypeError("str", s);
        return NULL;
    }
    Py_INCREF(s);
    PyUnicode_InternInPlace(&s);
    return s;
}

//////////////////// abs_longlong.proto ////////////////////
//@requires: ModuleSetupCode.c::IncludeStdlibH

static CYTHON_INLINE PY_LONG_LONG __Pyx_abs_longlong(PY_LONG_LONG x) {
#if defined (__cplusplus) && __cplusplus >= 201103L
    return std::abs(x);
#elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    return llabs(x);
#elif defined (_MSC_VER)
    // abs() is defined for long, but 64-bits type on MSVC is long long.
    // Use MS-specific _abs64() instead, which returns the original (negative) value for abs(-MAX-1)
    return _abs64(x);
#elif defined (__GNUC__)
    // gcc or clang on 64 bit windows.
    return __builtin_llabs(x);
#else
    if (sizeof(PY_LONG_LONG) <= sizeof(Py_ssize_t))
        return __Pyx_sst_abs(x);
    return (x<0) ? -x : x;
#endif
}


//////////////////// py_abs.proto ////////////////////

#if CYTHON_USE_PYLONG_INTERNALS
static PyObject *__Pyx_PyLong_AbsNeg(PyObject *num);/*proto*/

#define __Pyx_PyNumber_Absolute(x) \
    ((likely(PyLong_CheckExact(x))) ? \
         (likely(__Pyx_PyLong_IsNonNeg(x)) ? __Pyx_NewRef(x) : __Pyx_PyLong_AbsNeg(x)) : \
         PyNumber_Absolute(x))

#else
#define __Pyx_PyNumber_Absolute(x)  PyNumber_Absolute(x)
#endif

//////////////////// py_abs ////////////////////

#if CYTHON_USE_PYLONG_INTERNALS
static PyObject *__Pyx_PyLong_AbsNeg(PyObject *n) {
#if PY_VERSION_HEX >= 0x030C00A7
    if (likely(__Pyx_PyLong_IsCompact(n))) {
        return PyLong_FromSize_t(__Pyx_PyLong_CompactValueUnsigned(n));
    }
#else
    if (likely(Py_SIZE(n) == -1)) {
        // digits are unsigned
        return PyLong_FromUnsignedLong(__Pyx_PyLong_Digits(n)[0]);
    }
#endif
#if CYTHON_COMPILING_IN_CPYTHON
    {
        PyObject *copy = _PyLong_Copy((PyLongObject*)n);
        if (likely(copy)) {
            #if PY_VERSION_HEX >= 0x030C00A7
            // clear the sign bits to set the sign from SIGN_NEGATIVE (2) to positive (0)
            ((PyLongObject*)copy)->long_value.lv_tag ^= ((PyLongObject*)copy)->long_value.lv_tag & _PyLong_SIGN_MASK;
            #else
            // negate the size to swap the sign
            __Pyx_SET_SIZE(copy, -Py_SIZE(copy));
            #endif
        }
        return copy;
    }
#else
    return PyNumber_Negative(n);
#endif
}
#endif


//////////////////// pow2.proto ////////////////////

#define __Pyx_PyNumber_Power2(a, b) PyNumber_Power(a, b, Py_None)


//////////////////// divmod_int.proto //////////////////

static const {{RETURN_TYPE}} __Pyx_divmod_ERROR_VALUE_{{CFUNC_SUFFIX}} = {-1, -1};

static CYTHON_INLINE {{RETURN_TYPE}} __Pyx_divmod_{{CFUNC_SUFFIX}}({{TYPE}} a, {{TYPE}} b); /*proto*/


//////////////////// divmod_int //////////////////

static CYTHON_INLINE {{RETURN_TYPE}} __Pyx_divmod_{{CFUNC_SUFFIX}}({{TYPE}} a, {{TYPE}} b) {
    // Python and C/C++ use different algorithms in calculating quotients and remainders.
    // This results in different answers between Python and C/C++
    // when the dividend is negative and the divisor is positive and vice versa.
    {{TYPE}} q, r;
    if (unlikely(b == 0)) {
        {{if NOGIL}}PyGILState_STATE gilstate = PyGILState_Ensure();{{endif}}
        PyErr_SetString(PyExc_ZeroDivisionError, "division by zero");
        {{if NOGIL}}PyGILState_Release(gilstate);{{endif}}
        return __Pyx_divmod_ERROR_VALUE_{{CFUNC_SUFFIX}};
    } else if (a == 0) {
        q = 0;
        r = 0;
    } else if ((a < 0) != (b < 0)) {
        // see CMath.c :: DivInt and ModInt utility code
        q = a / b;
        r = a - q * b;
        {{TYPE}} adapt_python = ((r != 0) & ((r < 0) ^ (b < 0)));
        q -= adapt_python;
        r += adapt_python * b;
    }
    else {
        q = a / b;
        r = a % b;
    }

    {{RETURN_TYPE}} c_result = {q, r};
    return c_result;
}


//////////////////// divmod_float.proto //////////////////

static const {{RETURN_TYPE}} __Pyx_divmod_ERROR_VALUE_{{CFUNC_SUFFIX}} = {-1.0, -1.0};

static CYTHON_INLINE {{RETURN_TYPE}} __Pyx_divmod_{{CFUNC_SUFFIX}}({{TYPE}} a, {{TYPE}} b); /*proto*/


//////////////////// divmod_float //////////////////

static CYTHON_INLINE {{RETURN_TYPE}} __Pyx_divmod_{{CFUNC_SUFFIX}}({{TYPE}} a, {{TYPE}} b) {
    // Python and C/C++ use different algorithms in calculating quotients and remainders.
    // This results in different answers between Python and C/C++
    // when the dividend is negative and the divisor is positive and vice versa.

    // Adapted from CPython 3.14: floatobject.c / _float_div_mod()

    {{TYPE}} q, r, div;

    if (unlikely(b == 0.0)) {
        {{if NOGIL}}PyGILState_STATE gilstate = PyGILState_Ensure();{{endif}}
        PyErr_SetString(PyExc_ZeroDivisionError, "division by zero");
        {{if NOGIL}}PyGILState_Release(gilstate);{{endif}}
        return __Pyx_divmod_ERROR_VALUE_{{CFUNC_SUFFIX}};
    }

    r = fmod{{MATH_SUFFIX}}(a, b);
    // fmod is typically exact, so a-mod is *mathematically* an
    // exact multiple of b.  But this is fp arithmetic, and fp
    // a - mod is an approximation; the result is that div may
    // not be an exact integral value after the division, although
    // it will always be very close to one.
    div = (a - r) / b;
    if (r) {
        // ensure the remainder has the same sign as the denominator
        if ((b < 0) != (r < 0)) {
            r += b;
            div -= 1.0;
        }
    }
    else {
        // the remainder is zero, and in the presence of signed zeroes
        // fmod returns different results across platforms; ensure
        // it has the same sign as the denominator.
        r = copysign{{MATH_SUFFIX}}(0.0, b);
    }
    // snap quotient to nearest integral value
    if (div) {
        q = floor{{MATH_SUFFIX}}(div);
        if (div - q > 0.5) {
            q += 1.0;
        }
    }
    else {
        // div is zero - get the same sign as the true quotient
        q = copysign{{MATH_SUFFIX}}(0.0, a / b); /* zero w/ sign of a/b */
    }

    {{RETURN_TYPE}} c_result = {q, r};
    return c_result;
}


//////////////////// int_pyucs4.proto ////////////////////

static CYTHON_INLINE int __Pyx_int_from_UCS4(Py_UCS4 uchar);

//////////////////// int_pyucs4 ////////////////////

static int __Pyx_int_from_UCS4(Py_UCS4 uchar) {
    // Fast path for ascii digits
    if (likely(uchar >= (Py_UCS4)'0' && uchar <= (Py_UCS4)'9')) {
        return uchar - (Py_UCS4)'0';
    }
#if CYTHON_COMPILING_IN_LIMITED_API
    PyObject *u = PyUnicode_FromOrdinal(uchar);
    if (unlikely(!u)) return -1;
    PyObject *l = PyObject_CallFunctionObjArgs((PyObject*)(&PyLong_Type), u, NULL);
    Py_DECREF(u);
    if (unlikely(!l)) return -1;
#if __PYX_LIMITED_VERSION_HEX >= 0x030d0000
    int result = PyLong_AsInt(l);
#else
    // just don't handle overflow - it's very difficult to see how we'll get it from
    // a single digit.
    int result = (int)PyLong_AsLong(l);
#endif
    Py_DECREF(l);
    return result;
#else
    int digit = Py_UNICODE_TODECIMAL(uchar);
    if (unlikely(digit < 0)) {
        PyErr_Format(PyExc_ValueError,
            "invalid literal for int() with base 10: '%c'",
            (int) uchar);
        return -1;
    }
    return digit;
#endif
}


//////////////////// float_pyucs4.proto ////////////////////

static CYTHON_INLINE double __Pyx_double_from_UCS4(Py_UCS4 uchar);

//////////////////// float_pyucs4 ////////////////////

static double __Pyx_double_from_UCS4(Py_UCS4 uchar) {
    // fast path for "just an ascii digit"
    if (likely(uchar >= (Py_UCS4)'0' && uchar <= (Py_UCS4)'9')) {
        return uchar - (Py_UCS4)'0';
    }
#if CYTHON_COMPILING_IN_LIMITED_API
    PyObject *u = PyUnicode_FromOrdinal(uchar);
    if (unlikely(!u)) return -1.0;
    PyObject *f = PyFloat_FromString(u);
    Py_DECREF(u);
    if (unlikely(!f)) return -1.0;
    double result = PyFloat_AsDouble(f);
    Py_DECREF(f);
    return result;
#else
    // ...TONUMERIC would initially seem to be a better fit.
    // However, that accepts things like the "half" symbol, while
    // float(string) rejects those.
    double digit = Py_UNICODE_TODECIMAL(uchar);
    if (unlikely(digit < 0.0)) {
        PyErr_Format(PyExc_ValueError,
            "could not convert string to float: '%c'",
            (int) uchar);
        return -1.0;
    }
    return digit;
#endif
}


//////////////////// object_ord.proto ////////////////////
//@requires: TypeConversion.c::UnicodeAsUCS4

#define __Pyx_PyObject_Ord(c) \
    (likely(PyUnicode_Check(c)) ? (long)__Pyx_PyUnicode_AsPy_UCS4(c) : __Pyx__PyObject_Ord(c))
static long __Pyx__PyObject_Ord(PyObject* c); /*proto*/

//////////////////// object_ord ////////////////////

static long __Pyx__PyObject_Ord(PyObject* c) {
    Py_ssize_t size;
    if (PyBytes_Check(c)) {
        size = __Pyx_PyBytes_GET_SIZE(c);
        if (likely(size == 1)) {
#if CYTHON_ASSUME_SAFE_MACROS
            return (unsigned char) PyBytes_AS_STRING(c)[0];
#else
            char *data = PyBytes_AsString(c);
            if (unlikely(!data)) return -1;
            return (unsigned char) data[0];
#endif
        }
#if !CYTHON_ASSUME_SAFE_SIZE
        else if (unlikely(size < 0)) return -1;
#endif
    } else if (PyByteArray_Check(c)) {
        size = __Pyx_PyByteArray_GET_SIZE(c);
        if (likely(size == 1)) {
#if CYTHON_ASSUME_SAFE_MACROS
            return (unsigned char) PyByteArray_AS_STRING(c)[0];
#else
            char *data = PyByteArray_AsString(c);
            if (unlikely(!data)) return -1;
            return (unsigned char) data[0];
#endif
        }
#if !CYTHON_ASSUME_SAFE_SIZE
        else if (unlikely(size < 0)) return -1;
#endif
    } else {
        // FIXME: support character buffers - but CPython doesn't support them either
        __Pyx_TypeName c_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(c));
        PyErr_Format(PyExc_TypeError,
            "ord() expected string of length 1, but " __Pyx_FMT_TYPENAME " found",
            c_type_name);
        __Pyx_DECREF_TypeName(c_type_name);
        return (long)(Py_UCS4)-1;
    }
    PyErr_Format(PyExc_TypeError,
        "ord() expected a character, but string of length %zd found", size);
    return (long)(Py_UCS4)-1;
}


//////////////////// py_dict_keys.proto ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_Keys(PyObject* d); /*proto*/

//////////////////// py_dict_keys ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_Keys(PyObject* d) {
    return CALL_UNBOUND_METHOD(PyDict_Type, "keys", d);
}

//////////////////// py_dict_values.proto ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_Values(PyObject* d); /*proto*/

//////////////////// py_dict_values ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_Values(PyObject* d) {
    return CALL_UNBOUND_METHOD(PyDict_Type, "values", d);
}

//////////////////// py_dict_items.proto ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_Items(PyObject* d); /*proto*/

//////////////////// py_dict_items ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_Items(PyObject* d) {
    return CALL_UNBOUND_METHOD(PyDict_Type, "items", d);
}

//////////////////// py_dict_iterkeys.proto ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_IterKeys(PyObject* d); /*proto*/

//////////////////// py_dict_iterkeys ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_IterKeys(PyObject* d) {
    return CALL_UNBOUND_METHOD(PyDict_Type, "keys", d);
}

//////////////////// py_dict_itervalues.proto ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_IterValues(PyObject* d); /*proto*/

//////////////////// py_dict_itervalues ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_IterValues(PyObject* d) {
    return CALL_UNBOUND_METHOD(PyDict_Type, "values", d);
}

//////////////////// py_dict_iteritems.proto ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_IterItems(PyObject* d); /*proto*/

//////////////////// py_dict_iteritems ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_IterItems(PyObject* d) {
    return CALL_UNBOUND_METHOD(PyDict_Type, "items", d);
}

//////////////////// py_dict_viewkeys.proto ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_ViewKeys(PyObject* d); /*proto*/

//////////////////// py_dict_viewkeys ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_ViewKeys(PyObject* d) {
    return CALL_UNBOUND_METHOD(PyDict_Type, "keys", d);
}

//////////////////// py_dict_viewvalues.proto ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_ViewValues(PyObject* d); /*proto*/

//////////////////// py_dict_viewvalues ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_ViewValues(PyObject* d) {
    return CALL_UNBOUND_METHOD(PyDict_Type, "values", d);
}

//////////////////// py_dict_viewitems.proto ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_ViewItems(PyObject* d); /*proto*/

//////////////////// py_dict_viewitems ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyDict_ViewItems(PyObject* d) {
    return CALL_UNBOUND_METHOD(PyDict_Type, "items", d);
}


/////////////// dict_setdefault.proto ///////////////

static CYTHON_INLINE PyObject *__Pyx_PyDict_SetDefault(PyObject *d, PyObject *key, PyObject *default_value); /*proto*/

/////////////// dict_setdefault ///////////////

static CYTHON_INLINE PyObject *__Pyx_PyDict_SetDefault(PyObject *d, PyObject *key, PyObject *default_value) {
    PyObject* value;
#if __PYX_LIMITED_VERSION_HEX >= 0x030F0000 || (!CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX >= 0x030d00A4)
    PyDict_SetDefaultRef(d, key, default_value, &value);
#elif CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX >= 0x030C0000
    PyObject *args[] = {d, key, default_value};
    value = PyObject_VectorcallMethod(PYIDENT("setdefault"), args, 3 | PY_VECTORCALL_ARGUMENTS_OFFSET, NULL);
#elif CYTHON_COMPILING_IN_LIMITED_API
    value = PyObject_CallMethodObjArgs(d, PYIDENT("setdefault"), key, default_value, NULL);
#else
    value = PyDict_SetDefault(d, key, default_value);
    if (unlikely(!value)) return NULL;
    Py_INCREF(value);
#endif
    return value;
}


//////////////////// pyfrozenset_new.proto ////////////////////

static CYTHON_INLINE PyObject* __Pyx_PyFrozenSet_New(PyObject* it);

//////////////////// pyfrozenset_new ////////////////////
//@requires: ObjectHandling.c::PyObjectCallNoArg

static CYTHON_INLINE PyObject* __Pyx_PyFrozenSet_New(PyObject* it) {
    if (it) {
        PyObject* result;
#if CYTHON_COMPILING_IN_PYPY
        // PyPy currently lacks PyFrozenSet_CheckExact() and PyFrozenSet_New()
        PyObject* args;
        args = PyTuple_Pack(1, it);
        if (unlikely(!args))
            return NULL;
        result = PyObject_Call((PyObject*)&PyFrozenSet_Type, args, NULL);
        Py_DECREF(args);
        return result;
#else
        if (PyFrozenSet_CheckExact(it)) {
            Py_INCREF(it);
            return it;
        }
        result = PyFrozenSet_New(it);
        if (unlikely(!result))
            return NULL;
        if ((__PYX_LIMITED_VERSION_HEX >= 0x030A0000)
#if CYTHON_COMPILING_IN_LIMITED_API
            || __Pyx_get_runtime_version() >= 0x030A0000
#endif
            )
            return result;
        {
            Py_ssize_t size = __Pyx_PySet_GET_SIZE(result);
            if (likely(size > 0))
                return result;
#if !CYTHON_ASSUME_SAFE_SIZE
            if (unlikely(size < 0)) {
                Py_DECREF(result);
                return NULL;
            }
#endif
        }
        // empty frozenset is a singleton (on Python <3.10)
        // seems wasteful, but CPython does the same
        Py_DECREF(result);
#endif
    }
    return __Pyx_PyObject_CallNoArg((PyObject*) &PyFrozenSet_Type);
}


//////////////////// PySet_Update.proto ////////////////////

static CYTHON_INLINE int __Pyx_PySet_Update(PyObject* set, PyObject* it); /*proto*/

//////////////////// PySet_Update ////////////////////

static CYTHON_INLINE int __Pyx_PySet_Update(PyObject* set, PyObject* it) {
    PyObject *retval;
    #if CYTHON_USE_TYPE_SLOTS && !CYTHON_COMPILING_IN_PYPY
    if (PyAnySet_Check(it)) {
        Py_ssize_t size = __Pyx_PySet_GET_SIZE(it);
        #if !CYTHON_ASSUME_SAFE_SIZE
        if (unlikely(size < 0)) return -1;
        #endif
        if (size == 0)
            return 0;
        // fast and safe case: CPython will update our result set and return it
        retval = PySet_Type.tp_as_number->nb_inplace_or(set, it);
        if (likely(retval == set)) {
            Py_DECREF(retval);
            return 0;
        }
        if (unlikely(!retval))
            return -1;
        // unusual result, fall through to set.update() call below
        Py_DECREF(retval);
    }
    #endif
    retval = CALL_UNBOUND_METHOD(PySet_Type, "update", set, it);
    if (unlikely(!retval)) return -1;
    Py_DECREF(retval);
    return 0;
}

//////////////////// PyRange_Check.proto ////////////////////

#if CYTHON_COMPILING_IN_PYPY && !defined(PyRange_Check)
  #define PyRange_Check(obj)  __Pyx_TypeCheck((obj), &PyRange_Type)
#endif


///////////////// memoryview_get_from_buffer.proto ////////////////////

#if !CYTHON_COMPILING_IN_LIMITED_API
#define __Pyx_PyMemoryView_Get_{{name}}(o) PyMemoryView_GET_BUFFER(o)->{{name}}
#else
{{py:
out_types = dict(
    ndim='int', readonly='int',
    len='Py_ssize_t', itemsize='Py_ssize_t')
}} // can't get format like this unfortunately. It's unicode via getattr
{{py: out_type = out_types[name]}}
static {{out_type}} __Pyx_PyMemoryView_Get_{{name}}(PyObject *obj); /* proto */
#endif

////////////// memoryview_get_from_buffer /////////////////////////

#if !CYTHON_COMPILING_IN_LIMITED_API
#else
{{py:
out_types = dict(
    ndim='int', readonly='int',
    len='Py_ssize_t', itemsize='Py_ssize_t')
}}
{{py: out_type = out_types[name]}}
static {{out_type}} __Pyx_PyMemoryView_Get_{{name}}(PyObject *obj) {
    {{out_type}} result;
    PyObject *attr = PyObject_GetAttr(obj, PYIDENT("{{name}}"));
    if (!attr) {
        goto bad;
    }
{{if out_type == 'int'}}
    // I'm not worrying about overflow here because
    // ultimately it comes from a C struct that's an int
    result = PyLong_AsLong(attr);
{{elif out_type == 'Py_ssize_t'}}
    result = PyLong_AsSsize_t(attr);
{{endif}}
    Py_DECREF(attr);
    return result;

    bad:
    Py_XDECREF(attr);
    return -1;
}
#endif

////////////// PySliceAccessors.proto /////////////////////////

#if CYTHON_COMPILING_IN_LIMITED_API
#define __Pyx_PySlice_Start(o) PyObject_GetAttr(o, PYIDENT("start"))
#define __Pyx_PySlice_Stop(o) PyObject_GetAttr(o, PYIDENT("stop"))
#define __Pyx_PySlice_Step(o) PyObject_GetAttr(o, PYIDENT("step"))
#elif CYTHON_COMPILING_IN_GRAAL && defined(GRAALPY_VERSION_NUM) && GRAALPY_VERSION_NUM > 0x19000000
// Graal defines it's own accessor functions
#define __Pyx_PySlice_Start(o) GraalPySlice_Start(o)
#define __Pyx_PySlice_Stop(o) GraalPySlice_Stop(o)
#define __Pyx_PySlice_Step(o) GraalPySlice_Step(o)
#elif CYTHON_COMPILING_IN_GRAAL
// Remove when GraalPy 24 goes EOL
#define __Pyx_PySlice_Start(o) __Pyx_NewRef(PySlice_Start((PySliceObject*)o))
#define __Pyx_PySlice_Stop(o) __Pyx_NewRef(PySlice_Stop((PySliceObject*)o))
#define __Pyx_PySlice_Step(o) __Pyx_NewRef(PySlice_Step((PySliceObject*)o))
#else
#define __Pyx_PySlice_Start(o) __Pyx_NewRef(((PySliceObject*)o)->start)
#define __Pyx_PySlice_Stop(o) __Pyx_NewRef(((PySliceObject*)o)->stop)
#define __Pyx_PySlice_Step(o) __Pyx_NewRef(((PySliceObject*)o)->step)
#endif
