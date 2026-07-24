/*
These functions provide integer arithmetic with integer checking.  They do not
actually raise an exception when an overflow is detected, but rather set a bit
in the overflow parameter.  (This parameter may be reused across several
arithmetic operations, so should be or-ed rather than assigned to.)

The implementation is divided into two parts, the signed and unsigned basecases,
which is where the magic happens, and a generic template matching a specific
type to an implementation based on its (c-compile-time) size and signedness.

When possible, branching is avoided, and preference is given to speed over
accuracy (a low rate of falsely "detected" overflows are acceptable,
undetected overflows are not).


TODO: Hook up checking.
TODO: Conditionally support 128-bit with intmax_t?
*/

/////////////// Common.proto ///////////////

static int __Pyx_check_twos_complement(void) {
    if ((-1) != (~0)) {
        PyErr_SetString(PyExc_RuntimeError, "Two's complement required for overflow checks.");
        return 1;
    } else if ((sizeof(short) == sizeof(int))) {
        PyErr_SetString(PyExc_RuntimeError, "sizeof(short) < sizeof(int) required for overflow checks.");
        return 1;
    } else {
        return 0;
    }
}

#define __PYX_SIGN_BIT(type)    ((((unsigned type) 1) << (sizeof(type) * 8 - 1)))
#define __PYX_HALF_MAX(type)    ((((type) 1) << (sizeof(type) * 8 - 2)))
#define __PYX_MIN(type)         ((__PYX_IS_UNSIGNED(type) ? (type) 0 : 0 - __PYX_HALF_MAX(type) - __PYX_HALF_MAX(type)))
#define __PYX_MAX(type)         ((~__PYX_MIN(type)))

#define __Pyx_add_no_overflow(a, b, overflow) ((a) + (b))
#define __Pyx_add_const_no_overflow(a, b, overflow) ((a) + (b))
#define __Pyx_sub_no_overflow(a, b, overflow) ((a) - (b))
#define __Pyx_sub_const_no_overflow(a, b, overflow) ((a) - (b))
#define __Pyx_mul_no_overflow(a, b, overflow) ((a) * (b))
#define __Pyx_mul_const_no_overflow(a, b, overflow) ((a) * (b))
#define __Pyx_div_no_overflow(a, b, overflow) ((a) / (b))
#define __Pyx_div_const_no_overflow(a, b, overflow) ((a) / (b))

#if defined(__has_builtin)
#  if __has_builtin(__builtin_add_overflow) && !defined(__ibmxl__)
#    define __PYX_HAVE_BUILTIN_OVERFLOW
#  endif
#elif defined(__GNUC__) && (__GNUC__ >= 5) && (!defined(__INTEL_COMPILER) || (__INTEL_COMPILER >= 1800))
#  define __PYX_HAVE_BUILTIN_OVERFLOW
#endif

#if defined(__GNUC__)
#  define __Pyx_is_constant(x) (__builtin_constant_p(x))
#elif defined(__has_builtin)
#  if __has_builtin(__builtin_constant_p)
#    define __Pyx_is_constant(x) (__builtin_constant_p(x))
#  endif
#else
#  define __Pyx_is_constant(x) (0)
#endif

/////////////// Common.init ///////////////

if (likely(__Pyx_check_twos_complement() == 0)); else
// error propagation code is appended automatically

/////////////// BaseCaseUnsigned.proto ///////////////

static CYTHON_INLINE {{UINT}} __Pyx_add_{{NAME}}_checking_overflow({{UINT}} a, {{UINT}} b, int *overflow);
static CYTHON_INLINE {{UINT}} __Pyx_sub_{{NAME}}_checking_overflow({{UINT}} a, {{UINT}} b, int *overflow);
static CYTHON_INLINE {{UINT}} __Pyx_mul_{{NAME}}_checking_overflow({{UINT}} a, {{UINT}} b, int *overflow);
static CYTHON_INLINE {{UINT}} __Pyx_div_{{NAME}}_checking_overflow({{UINT}} a, {{UINT}} b, int *overflow);

// Use these when b is known at compile time.
#define __Pyx_add_const_{{NAME}}_checking_overflow __Pyx_add_{{NAME}}_checking_overflow
#define __Pyx_sub_const_{{NAME}}_checking_overflow __Pyx_sub_{{NAME}}_checking_overflow
#if defined(__PYX_HAVE_BUILTIN_OVERFLOW)
#define __Pyx_mul_const_{{NAME}}_checking_overflow __Pyx_mul_{{NAME}}_checking_overflow
#else
static CYTHON_INLINE {{UINT}} __Pyx_mul_const_{{NAME}}_checking_overflow({{UINT}} a, {{UINT}} constant, int *overflow);
#endif
#define __Pyx_div_const_{{NAME}}_checking_overflow __Pyx_div_{{NAME}}_checking_overflow

/////////////// BaseCaseUnsigned ///////////////

#if defined(__PYX_HAVE_BUILTIN_OVERFLOW)

static CYTHON_INLINE {{UINT}} __Pyx_add_{{NAME}}_checking_overflow({{UINT}} a, {{UINT}} b, int *overflow) {
    {{UINT}} result;
    *overflow |= __builtin_add_overflow(a, b, &result);
    return result;
}

static CYTHON_INLINE {{UINT}} __Pyx_sub_{{NAME}}_checking_overflow({{UINT}} a, {{UINT}} b, int *overflow) {
    {{UINT}} result;
    *overflow |= __builtin_sub_overflow(a, b, &result);
    return result;
}

static CYTHON_INLINE {{UINT}} __Pyx_mul_{{NAME}}_checking_overflow({{UINT}} a, {{UINT}} b, int *overflow) {
    {{UINT}} result;
    *overflow |= __builtin_mul_overflow(a, b, &result);
    return result;
}

#else

static CYTHON_INLINE {{UINT}} __Pyx_add_{{NAME}}_checking_overflow({{UINT}} a, {{UINT}} b, int *overflow) {
    {{UINT}} r = a + b;
    *overflow |= r < a;
    return r;
}

static CYTHON_INLINE {{UINT}} __Pyx_sub_{{NAME}}_checking_overflow({{UINT}} a, {{UINT}} b, int *overflow) {
    {{UINT}} r = a - b;
    *overflow |= r > a;
    return r;
}

static CYTHON_INLINE {{UINT}} __Pyx_mul_{{NAME}}_checking_overflow({{UINT}} a, {{UINT}} b, int *overflow) {
    // if we have a constant, use the constant version
    if (__Pyx_is_constant(b)) {
        return __Pyx_mul_const_{{NAME}}_checking_overflow(a, b, overflow);
    } else if (__Pyx_is_constant(a)) {
        return __Pyx_mul_const_{{NAME}}_checking_overflow(b, a, overflow);
    } else if ((sizeof({{UINT}}) < sizeof(unsigned long))) {
        unsigned long big_r = ((unsigned long) a) * ((unsigned long) b);
        {{UINT}} r = ({{UINT}}) big_r;
        *overflow |= big_r != r;
        return r;
    } else if ((sizeof({{UINT}}) < sizeof(unsigned PY_LONG_LONG))) {
        unsigned PY_LONG_LONG big_r = ((unsigned PY_LONG_LONG) a) * ((unsigned PY_LONG_LONG) b);
        {{UINT}} r = ({{UINT}}) big_r;
        *overflow |= big_r != r;
        return r;
    } else {
        return __Pyx_mul_const_{{NAME}}_checking_overflow(a, b, overflow);
    }
}

static CYTHON_INLINE {{UINT}} __Pyx_mul_const_{{NAME}}_checking_overflow({{UINT}} a, {{UINT}} b, int *overflow) {
    // note that deliberately the overflow check is written such that it divides by b; this
    // function is used when b is a constant thus the compiler should be able to eliminate the
    // (very slow on most CPUs!) division operation
    {{UINT}} prod;
    if (__Pyx_is_constant(a) && !__Pyx_is_constant(b)) {
        // if a is a compile-time constant and b isn't, swap them
        {{UINT}} temp = b;
        b = a;
        a = temp;
    }
    prod = a * b;
    if (b != 0)
        *overflow |= a > (__PYX_MAX({{UINT}}) / b);
    return prod;
}
#endif // __PYX_HAVE_BUILTIN_OVERFLOW


static CYTHON_INLINE {{UINT}} __Pyx_div_{{NAME}}_checking_overflow({{UINT}} a, {{UINT}} b, int *overflow) {
    if (b == 0) {
        *overflow |= 1;
        return 0;
    }
    return a / b;
}


/////////////// BaseCaseSigned.proto ///////////////

static CYTHON_INLINE {{INT}} __Pyx_add_{{NAME}}_checking_overflow({{INT}} a, {{INT}} b, int *overflow);
static CYTHON_INLINE {{INT}} __Pyx_sub_{{NAME}}_checking_overflow({{INT}} a, {{INT}} b, int *overflow);
static CYTHON_INLINE {{INT}} __Pyx_mul_{{NAME}}_checking_overflow({{INT}} a, {{INT}} b, int *overflow);
static CYTHON_INLINE {{INT}} __Pyx_div_{{NAME}}_checking_overflow({{INT}} a, {{INT}} b, int *overflow);


// Use when b is known at compile time.
#define __Pyx_add_const_{{NAME}}_checking_overflow __Pyx_add_{{NAME}}_checking_overflow
#define __Pyx_sub_const_{{NAME}}_checking_overflow __Pyx_sub_{{NAME}}_checking_overflow
#if defined(__PYX_HAVE_BUILTIN_OVERFLOW)
#define __Pyx_mul_const_{{NAME}}_checking_overflow __Pyx_mul_{{NAME}}_checking_overflow
#else
static CYTHON_INLINE {{INT}} __Pyx_mul_const_{{NAME}}_checking_overflow({{INT}} a, {{INT}} constant, int *overflow);
#endif
#define __Pyx_div_const_{{NAME}}_checking_overflow __Pyx_div_{{NAME}}_checking_overflow

/////////////// BaseCaseSigned ///////////////

#if defined(__PYX_HAVE_BUILTIN_OVERFLOW)

static CYTHON_INLINE {{INT}} __Pyx_add_{{NAME}}_checking_overflow({{INT}} a, {{INT}} b, int *overflow) {
    {{INT}} result;
    *overflow |= __builtin_add_overflow(a, b, &result);
    return result;
}

static CYTHON_INLINE {{INT}} __Pyx_sub_{{NAME}}_checking_overflow({{INT}} a, {{INT}} b, int *overflow) {
    {{INT}} result;
    *overflow |= __builtin_sub_overflow(a, b, &result);
    return result;
}

static CYTHON_INLINE {{INT}} __Pyx_mul_{{NAME}}_checking_overflow({{INT}} a, {{INT}} b, int *overflow) {
    {{INT}} result;
    *overflow |= __builtin_mul_overflow(a, b, &result);
    return result;
}

#else

static CYTHON_INLINE {{INT}} __Pyx_add_{{NAME}}_checking_overflow({{INT}} a, {{INT}} b, int *overflow) {
    if ((sizeof({{INT}}) < sizeof(long))) {
        long big_r = ((long) a) + ((long) b);
        {{INT}} r = ({{INT}}) big_r;
        *overflow |= big_r != r;
        return r;
    } else if ((sizeof({{INT}}) < sizeof(PY_LONG_LONG))) {
        PY_LONG_LONG big_r = ((PY_LONG_LONG) a) + ((PY_LONG_LONG) b);
        {{INT}} r = ({{INT}}) big_r;
        *overflow |= big_r != r;
        return r;
    } else {
        // Signed overflow undefined, but unsigned overflow is well defined. Casting is
        // implementation-defined, but we assume two's complement (see __Pyx_check_twos_complement
        // above), and arithmetic in two's-complement is the same as unsigned arithmetic.
        unsigned {{INT}} r = (unsigned {{INT}}) a + (unsigned {{INT}}) b;
        // Overflow happened if the operands have the same sign, but the result
        // has opposite sign.
        *overflow |= (((unsigned {{INT}})a ^ r) & ((unsigned {{INT}})b ^ r)) >> (8 * sizeof({{INT}}) - 1);
        return ({{INT}}) r;
    }
}

static CYTHON_INLINE {{INT}} __Pyx_sub_{{NAME}}_checking_overflow({{INT}} a, {{INT}} b, int *overflow) {
    // Compilers don't handle widening as well in the subtraction case, so don't bother
    unsigned {{INT}} r = (unsigned {{INT}}) a - (unsigned {{INT}}) b;
    // Overflow happened if the operands differing signs, and the result
    // has opposite sign to a.
    *overflow |= (((unsigned {{INT}})a ^ (unsigned {{INT}})b) & ((unsigned {{INT}})a ^ r)) >> (8 * sizeof({{INT}}) - 1);
    return ({{INT}}) r;
}

static CYTHON_INLINE {{INT}} __Pyx_mul_{{NAME}}_checking_overflow({{INT}} a, {{INT}} b, int *overflow) {
    // if we have a constant, use the constant version
    if (__Pyx_is_constant(b)) {
        return __Pyx_mul_const_{{NAME}}_checking_overflow(a, b, overflow);
    } else if (__Pyx_is_constant(a)) {
        return __Pyx_mul_const_{{NAME}}_checking_overflow(b, a, overflow);
    } else if ((sizeof({{INT}}) < sizeof(long))) {
        long big_r = ((long) a) * ((long) b);
        {{INT}} r = ({{INT}}) big_r;
        *overflow |= big_r != r;
        return ({{INT}}) r;
    } else if ((sizeof({{INT}}) < sizeof(PY_LONG_LONG))) {
        PY_LONG_LONG big_r = ((PY_LONG_LONG) a) * ((PY_LONG_LONG) b);
        {{INT}} r = ({{INT}}) big_r;
        *overflow |= big_r != r;
        return ({{INT}}) r;
    } else {
        return __Pyx_mul_const_{{NAME}}_checking_overflow(a, b, overflow);
    }
}

static CYTHON_INLINE {{INT}} __Pyx_mul_const_{{NAME}}_checking_overflow({{INT}} a, {{INT}} b, int *overflow) {
    // note that deliberately all these comparisons are written such that they divide by b; this
    // function is used when b is a constant thus the compiler should be able to eliminate the
    // (very slow on most CPUs!) division operations
    if (__Pyx_is_constant(a) && !__Pyx_is_constant(b)) {
        // if a is a compile-time constant and b isn't, swap them
        {{INT}} temp = b;
        b = a;
        a = temp;
    }
    if (b > 1) {
        *overflow |= a > __PYX_MAX({{INT}}) / b;
        *overflow |= a < __PYX_MIN({{INT}}) / b;
    } else if (b == -1) {
        *overflow |= a == __PYX_MIN({{INT}});
    } else if (b < -1) {
        *overflow |= a > __PYX_MIN({{INT}}) / b;
        *overflow |= a < __PYX_MAX({{INT}}) / b;
    }
    return ({{INT}}) (((unsigned {{INT}})a) * ((unsigned {{INT}}) b));
}
#endif  // defined(__PYX_HAVE_BUILTIN_OVERFLOW)

static CYTHON_INLINE {{INT}} __Pyx_div_{{NAME}}_checking_overflow({{INT}} a, {{INT}} b, int *overflow) {
    if (b == 0) {
        *overflow |= 1;
        return 0;
    }
    *overflow |= a == __PYX_MIN({{INT}}) && b == -1;
    return ({{INT}}) ((unsigned {{INT}}) a / (unsigned {{INT}}) b);
}


/////////////// SizeCheck.init ///////////////

if (likely(__Pyx_check_sane_{{NAME}}() == 0)); else
// error propagation code is appended automatically

/////////////// SizeCheck.proto ///////////////

static int __Pyx_check_sane_{{NAME}}(void) {
    if (((sizeof({{TYPE}}) <= sizeof(int)) ||
            (sizeof({{TYPE}}) == sizeof(PY_LONG_LONG)) ||
            (sizeof({{TYPE}}) == sizeof(long)))) {
        return 0;
    } else {
        PyErr_Format(PyExc_RuntimeError, \
            "Bad size for int type %.{{max(60, len(TYPE))}}s: %d", "{{TYPE}}", (int) sizeof({{TYPE}}));
        return 1;
    }
}


/////////////// Binop.proto ///////////////

static CYTHON_INLINE {{TYPE}} __Pyx_{{BINOP}}_{{NAME}}_checking_overflow({{TYPE}} a, {{TYPE}} b, int *overflow);

/////////////// Binop ///////////////

static CYTHON_INLINE {{TYPE}} __Pyx_{{BINOP}}_{{NAME}}_checking_overflow({{TYPE}} a, {{TYPE}} b, int *overflow) {
    if ((sizeof({{TYPE}}) < sizeof(int))) {
        return __Pyx_{{BINOP}}_no_overflow(a, b, overflow);
    } else if (__PYX_IS_UNSIGNED({{TYPE}})) {
        if ((sizeof({{TYPE}}) == sizeof(unsigned int))) {
            return ({{TYPE}}) __Pyx_{{BINOP}}_unsigned_int_checking_overflow(a, b, overflow);
        } else if ((sizeof({{TYPE}}) == sizeof(unsigned long))) {
            return ({{TYPE}}) __Pyx_{{BINOP}}_unsigned_long_checking_overflow(a, b, overflow);
        } else if ((sizeof({{TYPE}}) == sizeof(unsigned PY_LONG_LONG))) {
            return ({{TYPE}}) __Pyx_{{BINOP}}_unsigned_long_long_checking_overflow(a, b, overflow);
        } else {
            Py_FatalError(
                "__Pyx_{{BINOP}}_{{NAME}}_checking_overflow({{TYPE}}) executed an unexpected code path. "
                "Please report this as a bug in Cython"
            );
            return 0; /* handled elsewhere */
        }
    } else {
        if ((sizeof({{TYPE}}) == sizeof(int))) {
            return ({{TYPE}}) __Pyx_{{BINOP}}_int_checking_overflow(a, b, overflow);
        } else if ((sizeof({{TYPE}}) == sizeof(long))) {
            return ({{TYPE}}) __Pyx_{{BINOP}}_long_checking_overflow(a, b, overflow);
        } else if ((sizeof({{TYPE}}) == sizeof(PY_LONG_LONG))) {
            return ({{TYPE}}) __Pyx_{{BINOP}}_long_long_checking_overflow(a, b, overflow);
        } else {
            Py_FatalError(
                "__Pyx_{{BINOP}}_{{NAME}}_checking_overflow({{TYPE}}) executed an unexpected code path. "
                "Please report this as a bug in Cython"
            );
            return 0; /* handled elsewhere */
        }
    }
}

/////////////// LeftShift.proto ///////////////

static CYTHON_INLINE {{TYPE}} __Pyx_lshift_{{NAME}}_checking_overflow({{TYPE}} a, {{TYPE}} b, int *overflow) {
    int overflow_check =
#if {{SIGNED}}
        (a < 0) || (b < 0) ||
#endif
        // the following must be a _logical_ OR as the RHS is undefined if the LHS is true
        (b >= ({{TYPE}}) (8 * sizeof({{TYPE}}))) || (a > (__PYX_MAX({{TYPE}}) >> b));
    if (overflow_check) {
        *overflow |= 1;
        return 0;
    } else {
        return a << b;
    }
}
#define __Pyx_lshift_const_{{NAME}}_checking_overflow __Pyx_lshift_{{NAME}}_checking_overflow


/////////////// UnaryNegOverflows.proto ///////////////

// from intobject.c
#define __Pyx_UNARY_NEG_WOULD_OVERFLOW(x)    \
        (((x) < 0) & ((unsigned long)(x) == 0-(unsigned long)(x)))
