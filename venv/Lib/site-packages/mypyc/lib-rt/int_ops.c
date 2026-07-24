// Int primitive operations (tagged arbitrary-precision integers)
//
// These are registered in mypyc.primitives.int_ops.

#include <Python.h>
#include "CPy.h"

#ifdef _MSC_VER
#include <intrin.h>
#endif

#ifndef _WIN32
// On 64-bit Linux and macOS, ssize_t and long are both 64 bits, and
// PyLong_FromLong is faster than PyLong_FromSsize_t, so use the faster one
#define CPyLong_FromSsize_t PyLong_FromLong
#else
// On 64-bit Windows, ssize_t is 64 bits but long is 32 bits, so we
// can't use the above trick
#define CPyLong_FromSsize_t PyLong_FromSsize_t
#endif

#if defined(__GNUC__) || defined(__clang__)
#  if defined(__x86_64__) || defined(_M_X64) || defined(__aarch64__) || (defined(__SIZEOF_POINTER__) && __SIZEOF_POINTER__ == 8)
#    define CPY_CLZ(x) __builtin_clzll((unsigned long long)(x))
#    define CPY_BITS 64
#  else
#    define CPY_CLZ(x) __builtin_clz((unsigned int)(x))
#    define CPY_BITS 32
#  endif
#endif


CPyTagged CPyTagged_FromSsize_t(Py_ssize_t value) {
    // We use a Python object if the value shifted left by 1 is too
    // large for Py_ssize_t
    if (unlikely(CPyTagged_TooBig(value))) {
        PyObject *object = PyLong_FromSsize_t(value);
        return ((CPyTagged)object) | CPY_INT_TAG;
    } else {
        return value << 1;
    }
}

CPyTagged CPyTagged_FromVoidPtr(void *ptr) {
    if ((uintptr_t)ptr > PY_SSIZE_T_MAX) {
        PyObject *object = PyLong_FromVoidPtr(ptr);
        return ((CPyTagged)object) | CPY_INT_TAG;
    } else {
        return CPyTagged_FromSsize_t((Py_ssize_t)ptr);
    }
}

CPyTagged CPyTagged_FromInt64(int64_t value) {
    if (unlikely(CPyTagged_TooBigInt64(value))) {
        PyObject *object = PyLong_FromLongLong(value);
        return ((CPyTagged)object) | CPY_INT_TAG;
    } else {
        return value << 1;
    }
}

PyObject *CPyTagged_AsObject(CPyTagged x) {
    PyObject *value;
    if (unlikely(CPyTagged_CheckLong(x))) {
        value = CPyTagged_LongAsObject(x);
        Py_INCREF(value);
    } else {
        value = CPyLong_FromSsize_t(CPyTagged_ShortAsSsize_t(x));
        if (value == NULL) {
            CPyError_OutOfMemory();
        }
    }
    return value;
}

PyObject *CPyTagged_StealAsObject(CPyTagged x) {
    PyObject *value;
    if (unlikely(CPyTagged_CheckLong(x))) {
        value = CPyTagged_LongAsObject(x);
    } else {
        value = CPyLong_FromSsize_t(CPyTagged_ShortAsSsize_t(x));
        if (value == NULL) {
            CPyError_OutOfMemory();
        }
    }
    return value;
}

Py_ssize_t CPyTagged_AsSsize_t(CPyTagged x) {
    if (likely(CPyTagged_CheckShort(x))) {
        return CPyTagged_ShortAsSsize_t(x);
    } else {
        return PyLong_AsSsize_t(CPyTagged_LongAsObject(x));
    }
}

CPy_NOINLINE
void CPyTagged_IncRef(CPyTagged x) {
    if (unlikely(CPyTagged_CheckLong(x))) {
        Py_INCREF(CPyTagged_LongAsObject(x));
    }
}

CPy_NOINLINE
void CPyTagged_DecRef(CPyTagged x) {
    if (unlikely(CPyTagged_CheckLong(x))) {
        Py_DECREF(CPyTagged_LongAsObject(x));
    }
}

CPy_NOINLINE
void CPyTagged_XDecRef(CPyTagged x) {
    if (unlikely(CPyTagged_CheckLong(x))) {
        Py_XDECREF(CPyTagged_LongAsObject(x));
    }
}

// Tagged int negation slow path, where the result may be a long integer
CPyTagged CPyTagged_Negate_(CPyTagged num) {
    PyObject *num_obj = CPyTagged_AsObject(num);
    PyObject *result = PyNumber_Negative(num_obj);
    if (result == NULL) {
        CPyError_OutOfMemory();
    }
    Py_DECREF(num_obj);
    return CPyTagged_StealFromObject(result);
}

// Tagged int addition slow path, where the result may be a long integer
CPyTagged CPyTagged_Add_(CPyTagged left, CPyTagged right) {
    PyObject *left_obj = CPyTagged_AsObject(left);
    PyObject *right_obj = CPyTagged_AsObject(right);
    PyObject *result = PyNumber_Add(left_obj, right_obj);
    if (result == NULL) {
        CPyError_OutOfMemory();
    }
    Py_DECREF(left_obj);
    Py_DECREF(right_obj);
    return CPyTagged_StealFromObject(result);
}

// Tagged int subtraction slow path, where the result may be a long integer
CPyTagged CPyTagged_Subtract_(CPyTagged left, CPyTagged right) {
    PyObject *left_obj = CPyTagged_AsObject(left);
    PyObject *right_obj = CPyTagged_AsObject(right);
    PyObject *result = PyNumber_Subtract(left_obj, right_obj);
    if (result == NULL) {
        CPyError_OutOfMemory();
    }
    Py_DECREF(left_obj);
    Py_DECREF(right_obj);
    return CPyTagged_StealFromObject(result);
}

// Tagged int multiplication slow path, where the result may be a long integer
CPyTagged CPyTagged_Multiply_(CPyTagged left, CPyTagged right) {
    PyObject *left_obj = CPyTagged_AsObject(left);
    PyObject *right_obj = CPyTagged_AsObject(right);
    PyObject *result = PyNumber_Multiply(left_obj, right_obj);
    if (result == NULL) {
        CPyError_OutOfMemory();
    }
    Py_DECREF(left_obj);
    Py_DECREF(right_obj);
    return CPyTagged_StealFromObject(result);
}

// Tagged int // slow path, where the result may be a long integer (or raise)
CPyTagged CPyTagged_FloorDivide_(CPyTagged left, CPyTagged right) {
    PyObject *left_obj = CPyTagged_AsObject(left);
    PyObject *right_obj = CPyTagged_AsObject(right);
    PyObject *result = PyNumber_FloorDivide(left_obj, right_obj);
    Py_DECREF(left_obj);
    Py_DECREF(right_obj);
    // Handle exceptions honestly because it could be ZeroDivisionError
    if (result == NULL) {
        return CPY_INT_TAG;
    } else {
        return CPyTagged_StealFromObject(result);
    }
}

// Tagged int % slow path, where the result may be a long integer (or raise)
CPyTagged CPyTagged_Remainder_(CPyTagged left, CPyTagged right) {
    PyObject *left_obj = CPyTagged_AsObject(left);
    PyObject *right_obj = CPyTagged_AsObject(right);
    PyObject *result = PyNumber_Remainder(left_obj, right_obj);
    Py_DECREF(left_obj);
    Py_DECREF(right_obj);
    // Handle exceptions honestly because it could be ZeroDivisionError
    if (result == NULL) {
        return CPY_INT_TAG;
    } else {
        return CPyTagged_StealFromObject(result);
    }
}

bool CPyTagged_IsEq_(CPyTagged left, CPyTagged right) {
    if (CPyTagged_CheckShort(right)) {
        return false;
    } else {
        PyObject *left_obj = CPyTagged_AsObject(left);
        PyObject *right_obj = CPyTagged_AsObject(right);
        int result = PyObject_RichCompareBool(left_obj, right_obj, Py_EQ);
        Py_DECREF(left_obj);
        Py_DECREF(right_obj);
        if (result == -1) {
            CPyError_OutOfMemory();
        }
        return result;
    }
}

bool CPyTagged_IsLt_(CPyTagged left, CPyTagged right) {
    PyObject *left_obj = CPyTagged_AsObject(left);
    PyObject *right_obj = CPyTagged_AsObject(right);
    int result = PyObject_RichCompareBool(left_obj, right_obj, Py_LT);
    Py_DECREF(left_obj);
    Py_DECREF(right_obj);
    if (result == -1) {
        CPyError_OutOfMemory();
    }
    return result;
}

PyObject *CPyLong_FromStrWithBase(PyObject *o, CPyTagged base) {
    Py_ssize_t base_size_t = CPyTagged_AsSsize_t(base);
    return PyLong_FromUnicodeObject(o, base_size_t);
}

PyObject *CPyLong_FromStr(PyObject *o) {
    CPyTagged base = CPyTagged_FromSsize_t(10);
    return CPyLong_FromStrWithBase(o, base);
}

CPyTagged CPyTagged_FromFloat(double f) {
    if (f < ((double)CPY_TAGGED_MAX + 1.0) && f > (CPY_TAGGED_MIN - 1.0)) {
        return (Py_ssize_t)f << 1;
    }
    PyObject *o = PyLong_FromDouble(f);
    if (o == NULL)
        return CPY_INT_TAG;
    return CPyTagged_StealFromObject(o);
}

PyObject *CPyBool_Str(bool b) {
    return PyObject_Str(b ? Py_True : Py_False);
}

// Bitwise op '&', '|' or '^' using the generic (slow) API
static CPyTagged GenericBitwiseOp(CPyTagged a, CPyTagged b, char op) {
    PyObject *aobj = CPyTagged_AsObject(a);
    PyObject *bobj = CPyTagged_AsObject(b);
    PyObject *r;
    if (op == '&') {
        r = PyNumber_And(aobj, bobj);
    } else if (op == '|') {
        r = PyNumber_Or(aobj, bobj);
    } else {
        r = PyNumber_Xor(aobj, bobj);
    }
    if (unlikely(r == NULL)) {
        CPyError_OutOfMemory();
    }
    Py_DECREF(aobj);
    Py_DECREF(bobj);
    return CPyTagged_StealFromObject(r);
}

// Return pointer to digits of a PyLong object. If it's a short
// integer, place digits in the buffer buf instead to avoid memory
// allocation (it's assumed to be big enough). Return the number of
// digits in *size. *size is negative if the integer is negative.
static digit *GetIntDigits(CPyTagged n, Py_ssize_t *size, digit *buf) {
    if (CPyTagged_CheckShort(n)) {
        Py_ssize_t val = CPyTagged_ShortAsSsize_t(n);
        bool neg = val < 0;
        int len = 1;
        if (neg) {
            val = -val;
        }
        buf[0] = val & PyLong_MASK;
        if (val > (Py_ssize_t)PyLong_MASK) {
            val >>= PyLong_SHIFT;
            buf[1] = val & PyLong_MASK;
            if (val > (Py_ssize_t)PyLong_MASK) {
                buf[2] = val >> PyLong_SHIFT;
                len = 3;
            } else {
                len = 2;
            }
        }
        *size = neg ? -len : len;
        return buf;
    } else {
        PyLongObject *obj = (PyLongObject *)CPyTagged_LongAsObject(n);
        *size = CPY_LONG_SIZE_SIGNED(obj);
        return &CPY_LONG_DIGIT(obj, 0);
    }
}

// Shared implementation of bitwise '&', '|' and '^' (specified by op) for at least
// one long operand. This is somewhat optimized for performance.
CPyTagged CPyTagged_BitwiseLongOp_(CPyTagged a, CPyTagged b, char op) {
    // Directly access the digits, as there is no fast C API function for this.
    digit abuf[3];
    digit bbuf[3];
    Py_ssize_t asize;
    Py_ssize_t bsize;
    digit *adigits = GetIntDigits(a, &asize, abuf);
    digit *bdigits = GetIntDigits(b, &bsize, bbuf);

    if (unlikely(asize < 0 || bsize < 0)) {
        // Negative operand. This is slower, but bitwise ops on them are pretty rare.
        return GenericBitwiseOp(a, b, op);
    }
    // Optimized implementation for two non-negative integers.
    // Swap a and b as needed to ensure a is no longer than b.
    if (asize > bsize) {
        digit *tmp = adigits;
        adigits = bdigits;
        bdigits = tmp;
        Py_ssize_t tmp_size = asize;
        asize = bsize;
        bsize = tmp_size;
    }
    void *digits = NULL;
    PyLongWriter *writer = PyLongWriter_Create(0, op == '&' ? asize : bsize, &digits);
    if (unlikely(writer == NULL)) {
        CPyError_OutOfMemory();
    }
    Py_ssize_t i;
    if (op == '&') {
        for (i = 0; i < asize; i++) {
            ((digit *)digits)[i] = adigits[i] & bdigits[i];
        }
    } else {
        if (op == '|') {
            for (i = 0; i < asize; i++) {
                ((digit *)digits)[i] = adigits[i] | bdigits[i];
            }
        } else {
            for (i = 0; i < asize; i++) {
                ((digit *)digits)[i] = adigits[i] ^ bdigits[i];
            }
        }
        for (; i < bsize; i++) {
            ((digit *)digits)[i] = bdigits[i];
        }
    }
    return CPyTagged_StealFromObject(PyLongWriter_Finish(writer));
}

// Bitwise '~' slow path
CPyTagged CPyTagged_Invert_(CPyTagged num) {
    PyObject *obj = CPyTagged_AsObject(num);
    PyObject *result = PyNumber_Invert(obj);
    if (unlikely(result == NULL)) {
        CPyError_OutOfMemory();
    }
    Py_DECREF(obj);
    return CPyTagged_StealFromObject(result);
}

// Bitwise '>>' slow path
CPyTagged CPyTagged_Rshift_(CPyTagged left, CPyTagged right) {
    // Long integer or negative shift -- use generic op
    PyObject *lobj = CPyTagged_AsObject(left);
    PyObject *robj = CPyTagged_AsObject(right);
    PyObject *result = PyNumber_Rshift(lobj, robj);
    Py_DECREF(lobj);
    Py_DECREF(robj);
    if (result == NULL) {
        // Propagate error (could be negative shift count)
        return CPY_INT_TAG;
    }
    return CPyTagged_StealFromObject(result);
}

// Bitwise '<<' slow path
CPyTagged CPyTagged_Lshift_(CPyTagged left, CPyTagged right) {
    // Long integer or out of range shift -- use generic op
    PyObject *lobj = CPyTagged_AsObject(left);
    PyObject *robj = CPyTagged_AsObject(right);
    PyObject *result = PyNumber_Lshift(lobj, robj);
    Py_DECREF(lobj);
    Py_DECREF(robj);
    if (result == NULL) {
        // Propagate error (could be negative shift count)
        return CPY_INT_TAG;
    }
    return CPyTagged_StealFromObject(result);
}

// i64 unboxing slow path
int64_t CPyLong_AsInt64_(PyObject *o) {
    int overflow;
    int64_t result = PyLong_AsLongLongAndOverflow(o, &overflow);
    if (result == -1) {
        if (PyErr_Occurred()) {
            return CPY_LL_INT_ERROR;
        } else if (overflow) {
            PyErr_SetString(PyExc_ValueError, "int too large to convert to i64");
            return CPY_LL_INT_ERROR;
        }
    }
    return result;
}

int64_t CPyInt64_Divide(int64_t x, int64_t y) {
    if (y == 0) {
        PyErr_SetString(PyExc_ZeroDivisionError, "integer division or modulo by zero");
        return CPY_LL_INT_ERROR;
    }
    if (y == -1 && x == INT64_MIN) {
        PyErr_SetString(PyExc_OverflowError, "integer division overflow");
        return CPY_LL_INT_ERROR;
    }
    int64_t d = x / y;
    // Adjust for Python semantics
    if (((x < 0) != (y < 0)) && d * y != x) {
        d--;
    }
    return d;
}

int64_t CPyInt64_Remainder(int64_t x, int64_t y) {
    if (y == 0) {
        PyErr_SetString(PyExc_ZeroDivisionError, "integer division or modulo by zero");
        return CPY_LL_INT_ERROR;
    }
    // Edge case: avoid core dump
    if (y == -1 && x == INT64_MIN) {
        return 0;
    }
    int64_t d = x % y;
    // Adjust for Python semantics
    if (((x < 0) != (y < 0)) && d != 0) {
        d += y;
    }
    return d;
}

// i32 unboxing slow path
int32_t CPyLong_AsInt32_(PyObject *o) {
    int overflow;
    long result = PyLong_AsLongAndOverflow(o, &overflow);
    if (result > 0x7fffffffLL || result < -0x80000000LL) {
        overflow = 1;
        result = -1;
    }
    if (result == -1) {
        if (PyErr_Occurred()) {
            return CPY_LL_INT_ERROR;
        } else if (overflow) {
            PyErr_SetString(PyExc_ValueError, "int too large to convert to i32");
            return CPY_LL_INT_ERROR;
        }
    }
    return result;
}

int32_t CPyInt32_Divide(int32_t x, int32_t y) {
    if (y == 0) {
        PyErr_SetString(PyExc_ZeroDivisionError, "integer division or modulo by zero");
        return CPY_LL_INT_ERROR;
    }
    if (y == -1 && x == INT32_MIN) {
        PyErr_SetString(PyExc_OverflowError, "integer division overflow");
        return CPY_LL_INT_ERROR;
    }
    int32_t d = x / y;
    // Adjust for Python semantics
    if (((x < 0) != (y < 0)) && d * y != x) {
        d--;
    }
    return d;
}

int32_t CPyInt32_Remainder(int32_t x, int32_t y) {
    if (y == 0) {
        PyErr_SetString(PyExc_ZeroDivisionError, "integer division or modulo by zero");
        return CPY_LL_INT_ERROR;
    }
    // Edge case: avoid core dump
    if (y == -1 && x == INT32_MIN) {
        return 0;
    }
    int32_t d = x % y;
    // Adjust for Python semantics
    if (((x < 0) != (y < 0)) && d != 0) {
        d += y;
    }
    return d;
}

void CPyInt32_Overflow() {
    PyErr_SetString(PyExc_ValueError, "int too large to convert to i32");
}

// i16 unboxing slow path
int16_t CPyLong_AsInt16_(PyObject *o) {
    int overflow;
    long result = PyLong_AsLongAndOverflow(o, &overflow);
    if (result > 0x7fff || result < -0x8000) {
        overflow = 1;
        result = -1;
    }
    if (result == -1) {
        if (PyErr_Occurred()) {
            return CPY_LL_INT_ERROR;
        } else if (overflow) {
            PyErr_SetString(PyExc_ValueError, "int too large to convert to i16");
            return CPY_LL_INT_ERROR;
        }
    }
    return result;
}

int16_t CPyInt16_Divide(int16_t x, int16_t y) {
    if (y == 0) {
        PyErr_SetString(PyExc_ZeroDivisionError, "integer division or modulo by zero");
        return CPY_LL_INT_ERROR;
    }
    if (y == -1 && x == INT16_MIN) {
        PyErr_SetString(PyExc_OverflowError, "integer division overflow");
        return CPY_LL_INT_ERROR;
    }
    int16_t d = x / y;
    // Adjust for Python semantics
    if (((x < 0) != (y < 0)) && d * y != x) {
        d--;
    }
    return d;
}

int16_t CPyInt16_Remainder(int16_t x, int16_t y) {
    if (y == 0) {
        PyErr_SetString(PyExc_ZeroDivisionError, "integer division or modulo by zero");
        return CPY_LL_INT_ERROR;
    }
    // Edge case: avoid core dump
    if (y == -1 && x == INT16_MIN) {
        return 0;
    }
    int16_t d = x % y;
    // Adjust for Python semantics
    if (((x < 0) != (y < 0)) && d != 0) {
        d += y;
    }
    return d;
}

void CPyInt16_Overflow() {
    PyErr_SetString(PyExc_ValueError, "int too large to convert to i16");
}

// u8 unboxing slow path
uint8_t CPyLong_AsUInt8_(PyObject *o) {
    int overflow;
    long result = PyLong_AsLongAndOverflow(o, &overflow);
    if (result < 0 || result >= 256) {
        overflow = 1;
        result = -1;
    }
    if (result == -1) {
        if (PyErr_Occurred()) {
            return CPY_LL_UINT_ERROR;
        } else if (overflow) {
            PyErr_SetString(PyExc_ValueError, "int too large or small to convert to u8");
            return CPY_LL_UINT_ERROR;
        }
    }
    return result;
}

void CPyUInt8_Overflow() {
    PyErr_SetString(PyExc_ValueError, "int too large or small to convert to u8");
}

double CPyTagged_TrueDivide(CPyTagged x, CPyTagged y) {
    if (unlikely(y == 0)) {
        PyErr_SetString(PyExc_ZeroDivisionError, "division by zero");
        return CPY_FLOAT_ERROR;
    }
    if (likely(!CPyTagged_CheckLong(x) && !CPyTagged_CheckLong(y))) {
        return (double)((Py_ssize_t)x >> 1) / (double)((Py_ssize_t)y >> 1);
    } else {
        PyObject *xo = CPyTagged_AsObject(x);
        PyObject *yo = CPyTagged_AsObject(y);
        PyObject *result = PyNumber_TrueDivide(xo, yo);
        if (result == NULL) {
            Py_DECREF(xo);
            Py_DECREF(yo);
            return CPY_FLOAT_ERROR;
        }
        double value = PyFloat_AsDouble(result);
        Py_DECREF(xo);
        Py_DECREF(yo);
        Py_DECREF(result);
        return value;
    }
    return 1.0;
}

static PyObject *CPyLong_ToBytes(PyObject *v, Py_ssize_t length, int little_endian, int signed_flag) {
    // This is a wrapper for PyLong_AsByteArray and PyBytes_FromStringAndSize
    PyObject *result = PyBytes_FromStringAndSize(NULL, length);
    if (!result) {
        return NULL;
    }
    unsigned char *bytes = (unsigned char *)PyBytes_AS_STRING(result);
#if PY_VERSION_HEX >= 0x030D0000  // 3.13.0
    int res = _PyLong_AsByteArray((PyLongObject *)v, bytes, length, little_endian, signed_flag, 1);
#else
    int res = _PyLong_AsByteArray((PyLongObject *)v, bytes, length, little_endian, signed_flag);
#endif
    if (res < 0) {
        Py_DECREF(result);
        return NULL;
    }
    return result;
}

// int.to_bytes(length, byteorder, signed=False)
PyObject *CPyTagged_ToBytes(CPyTagged self, Py_ssize_t length, PyObject *byteorder, int signed_flag) {
    PyObject *pyint = CPyTagged_AsObject(self);
    if (!PyUnicode_Check(byteorder)) {
        Py_DECREF(pyint);
        PyErr_SetString(PyExc_TypeError, "byteorder must be str");
        return NULL;
    }
    const char *order = PyUnicode_AsUTF8(byteorder);
    if (!order) {
        Py_DECREF(pyint);
        return NULL;
    }
    int little_endian;
    if (strcmp(order, "big") == 0) {
        little_endian = 0;
    } else if (strcmp(order, "little") == 0) {
        little_endian = 1;
    } else {
        PyErr_SetString(PyExc_ValueError, "byteorder must be either 'little' or 'big'");
        return NULL;
    }
    PyObject *result = CPyLong_ToBytes(pyint, length, little_endian, signed_flag);
    Py_DECREF(pyint);
    return result;
}

// int.to_bytes(length, byteorder="little", signed=False)
PyObject *CPyTagged_ToLittleEndianBytes(CPyTagged self, Py_ssize_t length, int signed_flag) {
    PyObject *pyint = CPyTagged_AsObject(self);
    PyObject *result = CPyLong_ToBytes(pyint, length, 1, signed_flag);
    Py_DECREF(pyint);
    return result;
}

// int.to_bytes(length, "big", signed=False)
PyObject *CPyTagged_ToBigEndianBytes(CPyTagged self, Py_ssize_t length, int signed_flag) {
    PyObject *pyint = CPyTagged_AsObject(self);
    PyObject *result = CPyLong_ToBytes(pyint, length, 0, signed_flag);
    Py_DECREF(pyint);
    return result;
}

// int.bit_length()
CPyTagged CPyTagged_BitLength(CPyTagged self) {
    // Handle zero
    if (self == 0) {
        return 0;
    }

    // Fast path for small (tagged) ints
    if (CPyTagged_CheckShort(self)) {
        Py_ssize_t val = CPyTagged_ShortAsSsize_t(self);
        Py_ssize_t absval = val < 0 ? -val : val;
        int bits = 0;
        if (absval) {
#if defined(_MSC_VER)
    #if defined(_WIN64)
            unsigned long idx;
            if (_BitScanReverse64(&idx, (unsigned __int64)absval)) {
                bits = (int)(idx + 1);
            }
    #else
            unsigned long idx;
            if (_BitScanReverse(&idx, (unsigned long)absval)) {
                bits = (int)(idx + 1);
            }
    #endif
#elif defined(__GNUC__) || defined(__clang__)
            bits = (int)(CPY_BITS - CPY_CLZ(absval));
#else
            // Fallback to loop if no builtin
            while (absval) {
                absval >>= 1;
                bits++;
            }
#endif
        }
        return bits << 1;
    }

    // Slow path for big ints
    PyObject *pyint = CPyTagged_AsObject(self);
    int bits = _PyLong_NumBits(pyint);
    Py_DECREF(pyint);
    if (bits < 0) {
        // _PyLong_NumBits sets an error on failure
        return CPY_INT_TAG;
    }
    return bits << 1;
}
