
/////////////// CDivisionWarning.proto ///////////////

static int __Pyx_cdivision_warning(const char *, int); /* proto */

/////////////// CDivisionWarning ///////////////

static int __Pyx_cdivision_warning(const char *filename, int lineno) {
#if CYTHON_COMPILING_IN_PYPY
    // avoid compiler warnings
    filename++; lineno++;
    return PyErr_Warn(PyExc_RuntimeWarning,
                     "division with oppositely signed operands, C and Python semantics differ");
#else
    return PyErr_WarnExplicit(PyExc_RuntimeWarning,
                              "division with oppositely signed operands, C and Python semantics differ",
                              filename,
                              lineno,
                              __Pyx_MODULE_NAME,
                              NULL);
#endif
}


/////////////// DivInt.proto ///////////////

static CYTHON_INLINE %(type)s __Pyx_div_%(type_name)s(%(type)s, %(type)s, int b_is_constant); /* proto */

/////////////// DivInt ///////////////

static CYTHON_INLINE %(type)s __Pyx_div_%(type_name)s(%(type)s a, %(type)s b, int b_is_constant) {
    %(type)s q = a / b;
    %(type)s r = a - q*b;
    %(type)s adapt_python = (b_is_constant ?
        // Take advantage of constant folding for 'b'.
        ((r != 0) & ((r < 0) ^ (b < 0))) :
        ((r != 0) & ((r ^ b) < 0))
    );
    return q - adapt_python;
}


/////////////// ModInt.proto ///////////////

static CYTHON_INLINE %(type)s __Pyx_mod_%(type_name)s(%(type)s, %(type)s, int b_is_constant); /* proto */

/////////////// ModInt ///////////////

static CYTHON_INLINE %(type)s __Pyx_mod_%(type_name)s(%(type)s a, %(type)s b, int b_is_constant) {
    %(type)s r = a %% b;
    %(type)s adapt_python = (b_is_constant ?
        // Take advantage of constant folding for 'b'.
        ((r != 0) & ((r < 0) ^ (b < 0))) :
        ((r != 0) & ((r ^ b) < 0))
    );
    return r + adapt_python * b;
}


/////////////// ModFloat.proto ///////////////

static CYTHON_INLINE %(type)s __Pyx_mod_%(type_name)s(%(type)s, %(type)s, int b_is_constant); /* proto */

/////////////// ModFloat ///////////////

static CYTHON_INLINE %(type)s __Pyx_mod_%(type_name)s(%(type)s a, %(type)s b, int b_is_constant) {
    CYTHON_UNUSED_VAR(b_is_constant);
    %(type)s r = fmod%(math_h_modifier)s(a, b);
    r += ((r != 0) & ((r < 0) ^ (b < 0))) * b;
    return r;
}


/////////////// IntPow.proto ///////////////

static CYTHON_INLINE %(type)s %(func_name)s(%(type)s, %(type)s); /* proto */

/////////////// IntPow ///////////////

static CYTHON_INLINE %(type)s %(func_name)s(%(type)s b, %(type)s e) {
    %(type)s t = b;
    switch (e) {
        case 3:
            t *= b;
        CYTHON_FALLTHROUGH;
        case 2:
            t *= b;
        CYTHON_FALLTHROUGH;
        case 1:
            return t;
        case 0:
            return 1;
    }
    #if %(signed)s
    if (unlikely(e<0)) return 0;
    #endif
    t = 1;
    while (likely(e)) {
        t *= (b * (e&1)) | ((~e)&1);    /* 1 or b */
        b *= b;
        e >>= 1;
    }
    return t;
}
