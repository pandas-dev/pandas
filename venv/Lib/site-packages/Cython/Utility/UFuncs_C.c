///////////////////////// UFuncsInit.proto /////////////////////////
//@proto_block: utility_code_proto_before_types

#include <numpy/arrayobject.h>
#include <numpy/ufuncobject.h>

// account for change in type of arguments to PyUFuncGenericFunction in Numpy 1.19.x
// Unfortunately we can only test against Numpy version 1.20.x since it wasn't marked
// as an API break. Therefore, I'm "solving" the issue by casting function pointer types
// on lower Numpy versions.
#if NPY_API_VERSION >= 0x0000000e // Numpy 1.20.x
#define __PYX_PYUFUNCGENERICFUNCTION_CAST(x) x
#else
#define __PYX_PYUFUNCGENERICFUNCTION_CAST(x) (PyUFuncGenericFunction)x
#endif

/////////////////////// UFuncTypeHandling.proto //////////////

#define __PYX_GET_NPY_COMPLEX_TYPE(tp) \
    sizeof(tp) == sizeof(npy_cfloat) ? NPY_CFLOAT : \
    sizeof(tp) == sizeof(npy_cdouble) ? NPY_CDOUBLE : \
    sizeof(tp) == sizeof(npy_clongdouble) ? NPY_CLONGDOUBLE : NPY_NOTYPE

#define __PYX_GET_NPY_FLOAT_TYPE(tp) \
    sizeof(tp) == sizeof(npy_float) ? NPY_FLOAT : \
    sizeof(tp) == sizeof(npy_double) ? NPY_DOUBLE : \
    sizeof(tp) == sizeof(npy_longdouble) ? NPY_LONGDOUBLE : NPY_NOTYPE

#define __PYX_GET_NPY_UINT_TYPE(tp) \
    sizeof(tp) == 1 ? NPY_UINT8 : \
    sizeof(tp) == 2 ? NPY_UINT16 : \
    sizeof(tp) == 4 ? NPY_UINT32 : \
    sizeof(tp) == 8 ? NPY_UINT64 : NPY_NOTYPE

#define __PYX_GET_NPY_SINT_TYPE(tp) \
    sizeof(tp) == 1 ? NPY_INT8 : \
    sizeof(tp) == 2 ? NPY_INT16 : \
    sizeof(tp) == 4 ? NPY_INT32 : \
    sizeof(tp) == 8 ? NPY_INT64 : NPY_NOTYPE

#define __PYX_GET_NPY_INT_TYPE(tp) ( \
    (((tp)-1) > (tp)0) ? \
        (__PYX_GET_NPY_UINT_TYPE(tp)) : \
        (__PYX_GET_NPY_SINT_TYPE(tp)) )

/////////////////////// UFuncTypedef.proto ///////////////////
//@requires: UFuncTypeHandling

enum {
    /*
      Deduce what to tell Numpy about the extern typedef '{{type_cname}}'
    */
    __Pyx_typedef_ufunc_{{type_substituted_cname}} = {{macro_name}}({{type_cname}}),

    /* Validation:
       If the line below is failing then you are using the extern typedef
       '{{type_cname}}' as a parameter in a ufunc and Cython can't work out
       how to map the type to a Numpy type.
    */
    __Pyx_typedef_ufunc_validation_{{type_substituted_cname}} =
        1 / (int)(__Pyx_typedef_ufunc_{{type_substituted_cname}} - NPY_NOTYPE)
};

/////////////////////// UFuncConsts.proto ////////////////////

// getter functions because we can't forward-declare arrays
static PyUFuncGenericFunction* {{ufunc_funcs_name}}(void); /* proto */
static char* {{ufunc_types_name}}(void); /* proto */
static void* {{ufunc_data_name}}[] = {NULL};  /* always null */

/////////////////////// UFuncConsts /////////////////////////

static PyUFuncGenericFunction* {{ufunc_funcs_name}}(void) {
    static PyUFuncGenericFunction arr[] = {
        {{for loop, cname in looper(func_cnames)}}
        __PYX_PYUFUNCGENERICFUNCTION_CAST(&{{cname}}){{if not loop.last}},{{endif}}
        {{endfor}}
    };
    return arr;
}

static char* {{ufunc_types_name}}(void) {
    static char arr[] = {
        {{for loop, tp in looper(type_constants)}}
        {{tp}}{{if not loop.last}},{{endif}}
        {{endfor}}
    };
    return arr;
}
