#ifndef VEC_H_INCL
#define VEC_H_INCL

// Header for the implementation of librt.vecs, which defines the 'vec' type.
// Refer to librt_vecs.c for more detailed information.

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h>
#include <stdint.h>
#include "mypyc_util.h"

// ABI version -- only an exact match is compatible. This will only be changed in
// very exceptional cases (likely never) due to strict backward compatibility
// requirements.
#define LIBRT_VECS_ABI_VERSION 1

// API version -- more recent versions must maintain backward compatibility, i.e.
// we can add new features but not remove or change existing features (unless
// ABI version is changed, but see the comment above).
#define LIBRT_VECS_API_VERSION 1

// Magic (native) integer return value on exception. Caller must also
// use PyErr_Occurred() since this overlaps with valid integer values.
#define MYPYC_INT_ERROR -113

// Item type constants for supported packed/specialized item types; must be
// even but not a multiple of 4 (2 + 4 * n). Each of these has a corresponding
// distinct implementation C extension class. For example, vec[i64] has a
// different runtime type than vec[i32]. All other item types use generic
// implementations.
#define VEC_ITEM_TYPE_I64 2
#define VEC_ITEM_TYPE_I32 6
#define VEC_ITEM_TYPE_I16 10
#define VEC_ITEM_TYPE_U8  14
#define VEC_ITEM_TYPE_FLOAT 18
#define VEC_ITEM_TYPE_BOOL 22

static inline size_t Vec_IsMagicItemType(size_t item_type) {
    return item_type & 2;
}


// Buffer objects


// vecbuf[i64]
typedef struct _VecI64BufObject {
    PyObject_VAR_HEAD
    int64_t items[1];
} VecI64BufObject;

// vecbuf[i32]
typedef struct _VecI32BufObject {
    PyObject_VAR_HEAD
    int32_t items[1];
} VecI32BufObject;

// vecbuf[i16]
typedef struct _VecI16BufObject {
    PyObject_VAR_HEAD
    int16_t items[1];
} VecI16BufObject;

// vecbuf[u8]
typedef struct _VecU8BufObject {
    PyObject_VAR_HEAD
    uint8_t items[1];
} VecU8BufObject;

// vecbuf[float]
typedef struct _VecFloatBufObject {
    PyObject_VAR_HEAD
    double items[1];
} VecFloatBufObject;

// vecbuf[bool]
typedef struct _VecBoolBufObject {
    PyObject_VAR_HEAD
    char items[1];
} VecBoolBufObject;

// Simple generic vecbuf: vecbuf[t] when t is a type object
typedef struct _VecTBufObject {
    PyObject_VAR_HEAD
    // Tagged pointer to PyTypeObject *. The lowest bit is 1 for optional item type.
    size_t item_type;
    PyObject *items[1];
} VecTBufObject;

typedef struct _VecNestedBufItem {
    Py_ssize_t len;
    void *items;
} VecNestedBufItem;

// Nested vec type: vec[vec[...]], vec[vec[...] | None], etc.
typedef struct _VecNestedBufObject {
    PyObject_VAR_HEAD
    // Tagged pointer to PyTypeObject *. Lowest bit is set for optional item type.
    // The second lowest bit is set for a packed item type (VEC_ITEM_TYPE_*).
    size_t item_type;
    // Number of nested vec types (of any kind, at least 1)
    size_t depth;
    VecNestedBufItem items[1];
} VecNestedBufObject;


// Unboxed vec objects
//
// The items pointer points to the first element of the items array in the
// corresponding buffer object. Use VEC_*_BUF() to recover the buffer object
// from the items pointer (only needed on cold paths like grow/refcount).
// Items pointer is NULL for empty/uninitialized vecs.


typedef struct _VecI64 {
    Py_ssize_t len;
    int64_t *items;
} VecI64;

typedef struct _VecI32 {
    Py_ssize_t len;
    int32_t *items;
} VecI32;

typedef struct _VecI16 {
    Py_ssize_t len;
    int16_t *items;
} VecI16;

typedef struct _VecU8 {
    Py_ssize_t len;
    uint8_t *items;
} VecU8;

typedef struct _VecFloat {
    Py_ssize_t len;
    double *items;
} VecFloat;

typedef struct _VecBool {
    Py_ssize_t len;
    char *items;
} VecBool;

typedef struct _VecT {
    Py_ssize_t len;
    PyObject **items;
} VecT;

typedef struct _VecNested {
    Py_ssize_t len;
    VecNestedBufItem *items;
} VecNested;

// Recover buffer object from items pointer. Only valid when items != NULL.
#define VEC_I64_BUF_FROM_ITEMS(items_) \
    ((VecI64BufObject *)((char *)(items_) - offsetof(VecI64BufObject, items)))
#define VEC_I32_BUF_FROM_ITEMS(items_) \
    ((VecI32BufObject *)((char *)(items_) - offsetof(VecI32BufObject, items)))
#define VEC_I16_BUF_FROM_ITEMS(items_) \
    ((VecI16BufObject *)((char *)(items_) - offsetof(VecI16BufObject, items)))
#define VEC_U8_BUF_FROM_ITEMS(items_) \
    ((VecU8BufObject *)((char *)(items_) - offsetof(VecU8BufObject, items)))
#define VEC_FLOAT_BUF_FROM_ITEMS(items_) \
    ((VecFloatBufObject *)((char *)(items_) - offsetof(VecFloatBufObject, items)))
#define VEC_BOOL_BUF_FROM_ITEMS(items_) \
    ((VecBoolBufObject *)((char *)(items_) - offsetof(VecBoolBufObject, items)))
#define VEC_T_BUF_FROM_ITEMS(items_) \
    ((VecTBufObject *)((char *)(items_) - offsetof(VecTBufObject, items)))
#define VEC_NESTED_BUF_FROM_ITEMS(items_) \
    ((VecNestedBufObject *)((char *)(items_) - offsetof(VecNestedBufObject, items)))

#define VEC_I64_BUF(v) VEC_I64_BUF_FROM_ITEMS((v).items)
#define VEC_I32_BUF(v) VEC_I32_BUF_FROM_ITEMS((v).items)
#define VEC_I16_BUF(v) VEC_I16_BUF_FROM_ITEMS((v).items)
#define VEC_U8_BUF(v) VEC_U8_BUF_FROM_ITEMS((v).items)
#define VEC_FLOAT_BUF(v) VEC_FLOAT_BUF_FROM_ITEMS((v).items)
#define VEC_BOOL_BUF(v) VEC_BOOL_BUF_FROM_ITEMS((v).items)
#define VEC_T_BUF(v) VEC_T_BUF_FROM_ITEMS((v).items)
#define VEC_NESTED_BUF(v) VEC_NESTED_BUF_FROM_ITEMS((v).items)

// Type-specific incref/decref. Safe when items may be NULL.
#define VEC_I64_INCREF(v) do { if ((v).items) Py_INCREF(VEC_I64_BUF(v)); } while (0)
#define VEC_I64_DECREF(v) do { if ((v).items) Py_DECREF(VEC_I64_BUF(v)); } while (0)
#define VEC_I32_INCREF(v) do { if ((v).items) Py_INCREF(VEC_I32_BUF(v)); } while (0)
#define VEC_I32_DECREF(v) do { if ((v).items) Py_DECREF(VEC_I32_BUF(v)); } while (0)
#define VEC_I16_INCREF(v) do { if ((v).items) Py_INCREF(VEC_I16_BUF(v)); } while (0)
#define VEC_I16_DECREF(v) do { if ((v).items) Py_DECREF(VEC_I16_BUF(v)); } while (0)
#define VEC_U8_INCREF(v) do { if ((v).items) Py_INCREF(VEC_U8_BUF(v)); } while (0)
#define VEC_U8_DECREF(v) do { if ((v).items) Py_DECREF(VEC_U8_BUF(v)); } while (0)
#define VEC_FLOAT_INCREF(v) do { if ((v).items) Py_INCREF(VEC_FLOAT_BUF(v)); } while (0)
#define VEC_FLOAT_DECREF(v) do { if ((v).items) Py_DECREF(VEC_FLOAT_BUF(v)); } while (0)
#define VEC_BOOL_INCREF(v) do { if ((v).items) Py_INCREF(VEC_BOOL_BUF(v)); } while (0)
#define VEC_BOOL_DECREF(v) do { if ((v).items) Py_DECREF(VEC_BOOL_BUF(v)); } while (0)
#define VEC_T_INCREF(v) do { if ((v).items) Py_INCREF(VEC_T_BUF(v)); } while (0)
#define VEC_T_DECREF(v) do { if ((v).items) Py_DECREF(VEC_T_BUF(v)); } while (0)
#define VEC_NESTED_INCREF(v) do { if ((v).items) Py_INCREF(VEC_NESTED_BUF(v)); } while (0)
#define VEC_NESTED_DECREF(v) do { if ((v).items) Py_DECREF(VEC_NESTED_BUF(v)); } while (0)


// Boxed vec objects


// Arbitrary boxed vec object (only shared bits)
typedef struct _VecObject {
    PyObject_HEAD
    Py_ssize_t len;
} VecObject;

// Base vec type object (for isinstance checks)
// This is an abstract base type that all specialized vec types inherit from.
// It cannot be instantiated directly - only used for isinstance(x, vec).
typedef struct _VecBaseObject {
    PyObject_HEAD
} VecBaseObject;

// Boxed vec[i64]
typedef struct _VecI64Object {
    PyObject_HEAD
    VecI64 vec;
} VecI64Object;

// Boxed vec[i32]
typedef struct _VecI32Object {
    PyObject_HEAD
    VecI32 vec;
} VecI32Object;

// Boxed vec[i16]
typedef struct _VecI16Object {
    PyObject_HEAD
    VecI16 vec;
} VecI16Object;

// Boxed vec[u8]
typedef struct _VecU8Object {
    PyObject_HEAD
    VecU8 vec;
} VecU8Object;

// Boxed vec[float]
typedef struct _VecFloatObject {
    PyObject_HEAD
    VecFloat vec;
} VecFloatObject;

// Boxed vec[bool]
typedef struct _VecBoolObject {
    PyObject_HEAD
    VecBool vec;
} VecBoolObject;

// Simple boxed generic vecbuf: vecbuf[t] when t is a type object
typedef struct _VecTObject {
    PyObject_HEAD
    VecT vec;
} VecTObject;

// Extended generic vec type: vec[t | None], vec[vec[...]], etc.
typedef struct _VecNestedObject {
    PyObject_HEAD
    VecNested vec;
} VecNestedObject;


#ifndef MYPYC_DECLARED_tuple_T2V88
#define MYPYC_DECLARED_tuple_T2V88
typedef struct tuple_T2V88 {
    VecI64 f0;
    int64_t f1;
} tuple_T2V88;
static tuple_T2V88 tuple_undefined_T2V88 = { { -1, NULL } , 0 };
#endif

#ifndef MYPYC_DECLARED_tuple_T2V44
#define MYPYC_DECLARED_tuple_T2V44
typedef struct tuple_T2V44 {
    VecI32 f0;
    int32_t f1;
} tuple_T2V44;
static tuple_T2V44 tuple_undefined_T2V44 = { { -1, NULL } , 0 };
#endif

#ifndef MYPYC_DECLARED_tuple_T2V22
#define MYPYC_DECLARED_tuple_T2V22
typedef struct tuple_T2V22 {
    VecI16 f0;
    int16_t f1;
} tuple_T2V22;
static tuple_T2V22 tuple_undefined_T2V22 = { { -1, NULL } , 0 };
#endif

#ifndef MYPYC_DECLARED_tuple_T2VU1U1
#define MYPYC_DECLARED_tuple_T2VU1U1
typedef struct tuple_T2VU1U1 {
    VecU8 f0;
    uint8_t f1;
} tuple_T2VU1U1;
static tuple_T2VU1U1 tuple_undefined_T2VU1U1 = { { -1, NULL } , 0 };
#endif

#ifndef MYPYC_DECLARED_tuple_T2VFF
#define MYPYC_DECLARED_tuple_T2VFF
typedef struct tuple_T2VFF {
    VecFloat f0;
    double f1;
} tuple_T2VFF;
static tuple_T2VFF tuple_undefined_T2VFF = { { -1, NULL } , 0.0 };
#endif

#ifndef MYPYC_DECLARED_tuple_T2VCC
#define MYPYC_DECLARED_tuple_T2VCC
typedef struct tuple_T2VCC {
    VecBool f0;
    char f1;
} tuple_T2VCC;
static tuple_T2VCC tuple_undefined_T2VCC = { { -1, NULL } , 0 };
#endif

typedef tuple_T2V88 VecI64PopResult;
typedef tuple_T2V44 VecI32PopResult;
typedef tuple_T2V22 VecI16PopResult;
typedef tuple_T2VU1U1 VecU8PopResult;
typedef tuple_T2VFF VecFloatPopResult;
typedef tuple_T2VCC VecBoolPopResult;

// vec[i64] operations + type objects (stored in a capsule)
typedef struct _VecI64API {
    PyTypeObject *boxed_type;
    PyTypeObject *buf_type;
    VecI64 (*alloc)(Py_ssize_t, Py_ssize_t);
    PyObject *(*box)(VecI64);
    VecI64 (*unbox)(PyObject *);
    VecI64 (*convert_from_nested)(VecNestedBufItem);
    VecI64 (*append)(VecI64, int64_t);
    VecI64PopResult (*pop)(VecI64, Py_ssize_t);
    VecI64 (*remove)(VecI64, int64_t);
    // TODO: Py_ssize_t
    VecI64 (*slice)(VecI64, int64_t, int64_t);
    VecI64 (*from_iterable)(PyObject *, int64_t);
    VecI64 (*extend)(VecI64, PyObject *);
    VecI64 (*extend_vec)(VecI64, VecI64);
    PyObject *(*to_list)(VecI64);
    PyObject *(*to_tuple)(VecI64);
} VecI64API;

// vec[i32] operations + type objects (stored in a capsule)
typedef struct _VecI32API {
    PyTypeObject *boxed_type;
    PyTypeObject *buf_type;
    VecI32 (*alloc)(Py_ssize_t, Py_ssize_t);
    PyObject *(*box)(VecI32);
    VecI32 (*unbox)(PyObject *);
    VecI32 (*convert_from_nested)(VecNestedBufItem);
    VecI32 (*append)(VecI32, int32_t);
    VecI32PopResult (*pop)(VecI32, Py_ssize_t);
    VecI32 (*remove)(VecI32, int32_t);
    // TODO: Py_ssize_t
    VecI32 (*slice)(VecI32, int64_t, int64_t);
    VecI32 (*from_iterable)(PyObject *, int64_t);
    VecI32 (*extend)(VecI32, PyObject *);
    VecI32 (*extend_vec)(VecI32, VecI32);
    PyObject *(*to_list)(VecI32);
    PyObject *(*to_tuple)(VecI32);
} VecI32API;

// vec[i16] operations + type objects (stored in a capsule)
typedef struct _VecI16API {
    PyTypeObject *boxed_type;
    PyTypeObject *buf_type;
    VecI16 (*alloc)(Py_ssize_t, Py_ssize_t);
    PyObject *(*box)(VecI16);
    VecI16 (*unbox)(PyObject *);
    VecI16 (*convert_from_nested)(VecNestedBufItem);
    VecI16 (*append)(VecI16, int16_t);
    VecI16PopResult (*pop)(VecI16, Py_ssize_t);
    VecI16 (*remove)(VecI16, int16_t);
    // TODO: Py_ssize_t
    VecI16 (*slice)(VecI16, int64_t, int64_t);
    VecI16 (*from_iterable)(PyObject *, int64_t);
    VecI16 (*extend)(VecI16, PyObject *);
    VecI16 (*extend_vec)(VecI16, VecI16);
    PyObject *(*to_list)(VecI16);
    PyObject *(*to_tuple)(VecI16);
} VecI16API;

// vec[u8] operations + type objects (stored in a capsule)
typedef struct _VecU8API {
    PyTypeObject *boxed_type;
    PyTypeObject *buf_type;
    VecU8 (*alloc)(Py_ssize_t, Py_ssize_t);
    PyObject *(*box)(VecU8);
    VecU8 (*unbox)(PyObject *);
    VecU8 (*convert_from_nested)(VecNestedBufItem);
    VecU8 (*append)(VecU8, uint8_t);
    VecU8PopResult (*pop)(VecU8, Py_ssize_t);
    VecU8 (*remove)(VecU8, uint8_t);
    // TODO: Py_ssize_t
    VecU8 (*slice)(VecU8, int64_t, int64_t);
    VecU8 (*from_iterable)(PyObject *, int64_t);
    VecU8 (*extend)(VecU8, PyObject *);
    VecU8 (*extend_vec)(VecU8, VecU8);
    PyObject *(*to_list)(VecU8);
    PyObject *(*to_tuple)(VecU8);
} VecU8API;

// vec[float] operations + type objects (stored in a capsule)
typedef struct _VecFloatAPI {
    PyTypeObject *boxed_type;
    PyTypeObject *buf_type;
    VecFloat (*alloc)(Py_ssize_t, Py_ssize_t);
    PyObject *(*box)(VecFloat);
    VecFloat (*unbox)(PyObject *);
    VecFloat (*convert_from_nested)(VecNestedBufItem);
    VecFloat (*append)(VecFloat, double);
    VecFloatPopResult (*pop)(VecFloat, Py_ssize_t);
    VecFloat (*remove)(VecFloat, double);
    // TODO: Py_ssize_t
    VecFloat (*slice)(VecFloat, int64_t, int64_t);
    VecFloat (*from_iterable)(PyObject *, int64_t);
    VecFloat (*extend)(VecFloat, PyObject *);
    VecFloat (*extend_vec)(VecFloat, VecFloat);
    PyObject *(*to_list)(VecFloat);
    PyObject *(*to_tuple)(VecFloat);
} VecFloatAPI;

// vec[bool] operations + type objects (stored in a capsule)
typedef struct _VecBoolAPI {
    PyTypeObject *boxed_type;
    PyTypeObject *buf_type;
    VecBool (*alloc)(Py_ssize_t, Py_ssize_t);
    PyObject *(*box)(VecBool);
    VecBool (*unbox)(PyObject *);
    VecBool (*convert_from_nested)(VecNestedBufItem);
    VecBool (*append)(VecBool, char);
    VecBoolPopResult (*pop)(VecBool, Py_ssize_t);
    VecBool (*remove)(VecBool, char);
    // TODO: Py_ssize_t
    VecBool (*slice)(VecBool, int64_t, int64_t);
    VecBool (*from_iterable)(PyObject *, int64_t);
    VecBool (*extend)(VecBool, PyObject *);
    VecBool (*extend_vec)(VecBool, VecBool);
    PyObject *(*to_list)(VecBool);
    PyObject *(*to_tuple)(VecBool);
} VecBoolAPI;

#ifndef MYPYC_DECLARED_tuple_T2VOO
#define MYPYC_DECLARED_tuple_T2VOO
typedef struct tuple_T2VOO {
    VecT f0;
    PyObject *f1;
} tuple_T2VOO;
static tuple_T2VOO tuple_undefined_T2VOO = { { -1, NULL } , NULL };
#endif

typedef tuple_T2VOO VecTPopResult;

// vec[T] operations + type objects (stored in a capsule)
//
// T is a class type or class type | None
typedef struct _VecTAPI {
    PyTypeObject *boxed_type;
    PyTypeObject *buf_type;
    VecT (*alloc)(Py_ssize_t, Py_ssize_t, size_t);
    PyObject *(*box)(VecT, size_t);
    VecT (*unbox)(PyObject *, size_t);
    VecT (*convert_from_nested)(VecNestedBufItem);
    VecT (*append)(VecT, PyObject *, size_t);
    VecTPopResult (*pop)(VecT, Py_ssize_t);
    VecT (*remove)(VecT, PyObject *);
    // TODO: Py_ssize_t
    VecT (*slice)(VecT, int64_t, int64_t);
    VecT (*from_iterable)(size_t, PyObject *, int64_t);
    VecT (*extend)(VecT, PyObject *, size_t);
    VecT (*extend_vec)(VecT, VecT, size_t);
    PyObject *(*to_list)(VecT);
    PyObject *(*to_tuple)(VecT);
} VecTAPI;


#ifndef MYPYC_DECLARED_tuple_T2VvVi
#define MYPYC_DECLARED_tuple_T2VvVi
typedef struct tuple_T2VvVi {
    VecNested f0;
    VecNestedBufItem f1;
} tuple_T2VvVi;
static tuple_T2VvVi tuple_undefined_T2VvVi = { { -1, NULL } , { -1, 0 } };
#endif

typedef tuple_T2VvVi VecNestedPopResult;

// Nested vec operations + type objects (stored in a capsule)
typedef struct _VecNestedAPI {
    PyTypeObject *boxed_type;
    PyTypeObject *buf_type;
    VecNested (*alloc)(Py_ssize_t, Py_ssize_t, size_t, size_t depth);
    PyObject *(*box)(VecNested);
    VecNested (*unbox)(PyObject *, size_t, size_t depth);
    VecNested (*convert_from_nested)(VecNestedBufItem);
    VecNested (*append)(VecNested, VecNestedBufItem);
    VecNestedPopResult (*pop)(VecNested, Py_ssize_t);
    VecNested (*remove)(VecNested, VecNestedBufItem);
    // TODO: Py_ssize_t
    VecNested (*slice)(VecNested, int64_t, int64_t);
    VecNested (*extend)(VecNested, PyObject *);
    VecNested (*extend_vec)(VecNested, VecNested);
} VecNestedAPI;

typedef struct {
    int (*abi_version)(void);
    int (*api_version)(void);
    VecTAPI *t;
    VecNestedAPI *nested;
    VecI64API *i64;
    VecI32API *i32;
    VecI16API *i16;
    VecU8API *u8;
    VecFloatAPI *float_;
    VecBoolAPI *bool_;
    PyTypeObject *(*get_vec_type)(void);  // Function to get base VecType for isinstance checks
} VecCapsule;

#define VEC_BUF_SIZE(b) ((b)->ob_base.ob_size)
#define VEC_ITEM_TYPE(t) ((PyTypeObject *)((t) & ~1))
#define VEC_BUF_ITEM_TYPE(b) VEC_ITEM_TYPE((b)->item_type)
#define VEC_IS_ERROR(v) ((v).len < 0)

// Type objects

// Buffer type objects that store vec items
extern PyTypeObject VecI64BufType;
extern PyTypeObject VecI32BufType;
extern PyTypeObject VecI16BufType;
extern PyTypeObject VecU8BufType;
extern PyTypeObject VecFloatBufType;
extern PyTypeObject VecBoolBufType;
extern PyTypeObject VecTBufType;
extern PyTypeObject VecNestedBufType;

// Wrapper type objects for boxed vec values
extern PyTypeObject VecI64Type;
extern PyTypeObject VecI32Type;
extern PyTypeObject VecI16Type;
extern PyTypeObject VecU8Type;
extern PyTypeObject VecFloatType;
extern PyTypeObject VecBoolType;
extern PyTypeObject VecTType;
extern PyTypeObject VecNestedType;

// Iterator type objects for vec iteration
extern PyTypeObject VecI64IterType;
extern PyTypeObject VecI32IterType;
extern PyTypeObject VecI16IterType;
extern PyTypeObject VecU8IterType;
extern PyTypeObject VecFloatIterType;
extern PyTypeObject VecBoolIterType;
extern PyTypeObject VecTIterType;
extern PyTypeObject VecNestedIterType;

// Type objects corresponding to the 'i64', 'i32', 'i16, and 'u8' types
extern PyTypeObject *LibRTVecs_I64TypeObj;
extern PyTypeObject *LibRTVecs_I32TypeObj;
extern PyTypeObject *LibRTVecs_I16TypeObj;
extern PyTypeObject *LibRTVecs_U8TypeObj;

extern VecI64API Vec_I64API;
extern VecI32API Vec_I32API;
extern VecI16API Vec_I16API;
extern VecU8API Vec_U8API;
extern VecFloatAPI Vec_FloatAPI;
extern VecBoolAPI Vec_BoolAPI;
extern VecTAPI Vec_TAPI;
extern VecNestedAPI Vec_NestedAPI;

static inline int Vec_CheckFloatError(PyObject *o) {
    if (PyFloat_Check(o)) {
        PyErr_SetString(PyExc_TypeError, "integer argument expected, got float");
        return 1;
    }
    return 0;
}

// vec[i64] operations

static inline int VecI64_Check(PyObject *o) {
    return o->ob_type == &VecI64Type;
}

static inline PyObject *VecI64_BoxItem(int64_t x) {
    return PyLong_FromLongLong(x);
}

static inline int64_t VecI64_UnboxItem(PyObject *o) {
    if (Vec_CheckFloatError(o))
        return -1;
    return PyLong_AsLongLong(o);
}

static inline int VecI64_IsUnboxError(int64_t x) {
    return x == -1 && PyErr_Occurred();
}

PyObject *VecI64_Box(VecI64);
VecI64 VecI64_Append(VecI64, int64_t x);
VecI64 VecI64_Extend(VecI64, PyObject *iterable);
VecI64 VecI64_ExtendVec(VecI64 dst, VecI64 src);
VecI64 VecI64_Remove(VecI64, int64_t x);
VecI64PopResult VecI64_Pop(VecI64 v, Py_ssize_t index);

// vec[i32] operations

static inline int VecI32_Check(PyObject *o) {
    return o->ob_type == &VecI32Type;
}

static inline PyObject *VecI32_BoxItem(int32_t x) {
    return PyLong_FromLongLong(x);
}

static inline int32_t VecI32_UnboxItem(PyObject *o) {
    if (Vec_CheckFloatError(o))
        return -1;
    long x = PyLong_AsLong(o);
    if (x > INT32_MAX || x < INT32_MIN) {
        PyErr_SetString(PyExc_OverflowError, "Python int too large to convert to i32");
        return -1;
    }
    return x;
}

static inline int VecI32_IsUnboxError(int32_t x) {
    return x == -1 && PyErr_Occurred();
}

PyObject *VecI32_Box(VecI32);
VecI32 VecI32_Append(VecI32, int32_t x);
VecI32 VecI32_Extend(VecI32, PyObject *iterable);
VecI32 VecI32_ExtendVec(VecI32 dst, VecI32 src);
VecI32 VecI32_Remove(VecI32, int32_t x);
VecI32PopResult VecI32_Pop(VecI32 v, Py_ssize_t index);

// vec[i16] operations

static inline int VecI16_Check(PyObject *o) {
    return o->ob_type == &VecI16Type;
}

static inline PyObject *VecI16_BoxItem(int16_t x) {
    return PyLong_FromLongLong(x);
}

static inline int16_t VecI16_UnboxItem(PyObject *o) {
    if (Vec_CheckFloatError(o))
        return -1;
    long x = PyLong_AsLong(o);
    if (x >= 32768 || x < -32768) {
        PyErr_SetString(PyExc_OverflowError, "Python int too large to convert to i16");
        return -1;
    }
    return x;
}

static inline int VecI16_IsUnboxError(int16_t x) {
    return x == -1 && PyErr_Occurred();
}

PyObject *VecI16_Box(VecI16);
VecI16 VecI16_Append(VecI16, int16_t x);
VecI16 VecI16_Extend(VecI16, PyObject *iterable);
VecI16 VecI16_ExtendVec(VecI16 dst, VecI16 src);
VecI16 VecI16_Remove(VecI16, int16_t x);
VecI16PopResult VecI16_Pop(VecI16 v, Py_ssize_t index);

// vec[u8] operations

static inline int VecU8_Check(PyObject *o) {
    return o->ob_type == &VecU8Type;
}

static inline PyObject *VecU8_BoxItem(uint8_t x) {
    return PyLong_FromUnsignedLong(x);
}

static inline uint8_t VecU8_UnboxItem(PyObject *o) {
    if (Vec_CheckFloatError(o))
        return -1;
    unsigned long x = PyLong_AsUnsignedLong(o);
    if (x <= 255)
        return x;
    else if (x == (unsigned long)-1)
        return 239;
    else {
        PyErr_SetString(PyExc_OverflowError, "Python int too large to convert to u8");
        return 239;
    }
}

static inline int VecU8_IsUnboxError(uint8_t x) {
    return x == 239 && PyErr_Occurred();
}

PyObject *VecU8_Box(VecU8);
VecU8 VecU8_Append(VecU8, uint8_t x);
VecU8 VecU8_Extend(VecU8, PyObject *iterable);
VecU8 VecU8_ExtendVec(VecU8 dst, VecU8 src);
VecU8 VecU8_Remove(VecU8, uint8_t x);
VecU8PopResult VecU8_Pop(VecU8 v, Py_ssize_t index);

// vec[float] operations

static inline int VecFloat_Check(PyObject *o) {
    return o->ob_type == &VecFloatType;
}

static inline PyObject *VecFloat_BoxItem(double x) {
    return PyFloat_FromDouble(x);
}

static inline double VecFloat_UnboxItem(PyObject *o) {
    return PyFloat_AsDouble(o);
}

static inline int VecFloat_IsUnboxError(double x) {
    return x == -1.0 && PyErr_Occurred();
}

PyObject *VecFloat_Box(VecFloat);
VecFloat VecFloat_Append(VecFloat, double x);
VecFloat VecFloat_Extend(VecFloat, PyObject *iterable);
VecFloat VecFloat_ExtendVec(VecFloat dst, VecFloat src);
VecFloat VecFloat_Remove(VecFloat, double x);
VecFloatPopResult VecFloat_Pop(VecFloat v, Py_ssize_t index);

// vec[bool] operations

static inline int VecBool_Check(PyObject *o) {
    return o->ob_type == &VecBoolType;
}

static inline PyObject *VecBool_BoxItem(char x) {
    if (x == 1) {
        Py_INCREF(Py_True);
        return Py_True;
    } else {
        Py_INCREF(Py_False);
        return Py_False;
    }
}

static inline char VecBool_UnboxItem(PyObject *o) {
    if (o == Py_False) {
        return 0;
    } else if (o == Py_True) {
        return 1;
    } else {
        PyErr_SetString(PyExc_TypeError, "bool value expected");
        return 2;
    }
}

static inline int VecBool_IsUnboxError(char x) {
    return x == 2;
}

PyObject *VecBool_Box(VecBool);
VecBool VecBool_Append(VecBool, char x);
VecBool VecBool_Extend(VecBool, PyObject *iterable);
VecBool VecBool_ExtendVec(VecBool dst, VecBool src);
VecBool VecBool_Remove(VecBool, char x);
VecBoolPopResult VecBool_Pop(VecBool v, Py_ssize_t index);

// vec[t] operations

static inline int VecT_Check(PyObject *o) {
    return o->ob_type == &VecTType;
}

static inline int VecT_ItemCheck(VecT v, PyObject *item, size_t item_type) {
    if (PyObject_TypeCheck(item, VEC_ITEM_TYPE(item_type))) {
        return 1;
    } else if ((item_type & 1) && item == Py_None) {
        return 1;
    } else {
        // TODO: better error message
        PyErr_SetString(PyExc_TypeError, "invalid item type");
        return 0;
    }
}

VecT VecT_New(Py_ssize_t size, Py_ssize_t cap, size_t item_type);
VecT VecT_FromIterable(size_t item_type, PyObject *iterable, int64_t cap);
PyObject *VecT_Box(VecT vec, size_t item_type);
VecT VecT_Append(VecT vec, PyObject *x, size_t item_type);
VecT VecT_Extend(VecT vec, PyObject *iterable, size_t item_type);
VecT VecT_ExtendVec(VecT dst, VecT src, size_t item_type);
VecT VecT_Remove(VecT vec, PyObject *x);
VecTPopResult VecT_Pop(VecT v, Py_ssize_t index);

// Nested vec operations

static inline int VecNested_Check(PyObject *o) {
    return o->ob_type == &VecNestedType;
}

VecNested VecNested_New(Py_ssize_t size, Py_ssize_t cap, size_t item_type, size_t depth);
PyObject *VecNested_FromIterable(size_t item_type, size_t depth, PyObject *iterable, int64_t cap);
PyObject *VecNested_Box(VecNested);
VecNested VecNested_Append(VecNested vec, VecNestedBufItem x);
VecNested VecNested_Extend(VecNested vec, PyObject *iterable);
VecNested VecNested_ExtendVec(VecNested dst, VecNested src);
VecNested VecNested_Remove(VecNested vec, VecNestedBufItem x);
VecNestedPopResult VecNested_Pop(VecNested v, Py_ssize_t index);

static inline PyObject *VecNested_ItemBuf(VecNestedBufObject *parent, VecNestedBufItem item) {
    if (item.items == NULL)
        return NULL;
    if (parent->depth > 1)
        return (PyObject *)VEC_NESTED_BUF_FROM_ITEMS(item.items);
    if (parent->item_type == VEC_ITEM_TYPE_I64)
        return (PyObject *)VEC_I64_BUF_FROM_ITEMS(item.items);
    if (parent->item_type == VEC_ITEM_TYPE_U8)
        return (PyObject *)VEC_U8_BUF_FROM_ITEMS(item.items);
    if (parent->item_type == VEC_ITEM_TYPE_FLOAT)
        return (PyObject *)VEC_FLOAT_BUF_FROM_ITEMS(item.items);
    if (parent->item_type == VEC_ITEM_TYPE_I32)
        return (PyObject *)VEC_I32_BUF_FROM_ITEMS(item.items);
    if (parent->item_type == VEC_ITEM_TYPE_I16)
        return (PyObject *)VEC_I16_BUF_FROM_ITEMS(item.items);
    if (parent->item_type == VEC_ITEM_TYPE_BOOL)
        return (PyObject *)VEC_BOOL_BUF_FROM_ITEMS(item.items);
    return (PyObject *)VEC_T_BUF_FROM_ITEMS(item.items);
}

static inline void VecNested_ItemXINCREF(VecNestedBufObject *parent, VecNestedBufItem item) {
    Py_XINCREF(VecNested_ItemBuf(parent, item));
}

static inline void VecNested_ItemXDECREF(VecNestedBufObject *parent, VecNestedBufItem item) {
    Py_XDECREF(VecNested_ItemBuf(parent, item));
}

static inline void VecNested_ItemCLEAR(VecNestedBufObject *parent, VecNestedBufItem *item) {
    PyObject *buf = VecNested_ItemBuf(parent, *item);
    item->items = NULL;
    Py_XDECREF(buf);
}

static inline int VecNested_ItemVISIT(
        VecNestedBufObject *parent, VecNestedBufItem item, visitproc visit, void *arg) {
    PyObject *buf = VecNested_ItemBuf(parent, item);
    if (buf) {
        int ret = visit(buf, arg);
        if (ret)
            return ret;
    }
    return 0;
}

// Return 0 on success, -1 on error. Store unboxed item in *unboxed if successful.
// Return a *borrowed* reference.
static inline int VecNested_UnboxItem(VecNested v, PyObject *item, VecNestedBufItem *unboxed) {
    VecNestedBufObject *v_buf = VEC_NESTED_BUF(v);
    size_t depth = v_buf->depth;
    if (depth == 1) {
        if (item->ob_type == &VecTType) {
            // Boxed vec[t] always has items != NULL (buf is allocated on boxing)
            VecTObject *o = (VecTObject *)item;
            if (VEC_T_BUF(o->vec)->item_type == v_buf->item_type) {
                unboxed->len = o->vec.len;
                unboxed->items = o->vec.items;
                return 0;
            }
        } else if (item->ob_type == &VecI64Type && v_buf->item_type == VEC_ITEM_TYPE_I64) {
            VecI64Object *o = (VecI64Object *)item;
            unboxed->len = o->vec.len;
            unboxed->items = o->vec.items;
            return 0;
        } else if (item->ob_type == &VecU8Type && v_buf->item_type == VEC_ITEM_TYPE_U8) {
            VecU8Object *o = (VecU8Object *)item;
            unboxed->len = o->vec.len;
            unboxed->items = o->vec.items;
            return 0;
        } else if (item->ob_type == &VecFloatType && v_buf->item_type == VEC_ITEM_TYPE_FLOAT) {
            VecFloatObject *o = (VecFloatObject *)item;
            unboxed->len = o->vec.len;
            unboxed->items = o->vec.items;
            return 0;
        } else if (item->ob_type == &VecI32Type && v_buf->item_type == VEC_ITEM_TYPE_I32) {
            VecI32Object *o = (VecI32Object *)item;
            unboxed->len = o->vec.len;
            unboxed->items = o->vec.items;
            return 0;
        } else if (item->ob_type == &VecI16Type && v_buf->item_type == VEC_ITEM_TYPE_I16) {
            VecI16Object *o = (VecI16Object *)item;
            unboxed->len = o->vec.len;
            unboxed->items = o->vec.items;
            return 0;
        } else if (item->ob_type == &VecBoolType && v_buf->item_type == VEC_ITEM_TYPE_BOOL) {
            VecBoolObject *o = (VecBoolObject *)item;
            unboxed->len = o->vec.len;
            unboxed->items = o->vec.items;
            return 0;
        }
    } else if (item->ob_type == &VecNestedType) {
        VecNestedObject *o = (VecNestedObject *)item;
        VecNestedBufObject *o_buf = VEC_NESTED_BUF(o->vec);
        if (o_buf->depth == v_buf->depth - 1
            && o_buf->item_type == v_buf->item_type) {
            unboxed->len = o->vec.len;
            unboxed->items = o->vec.items;
            return 0;
        }
    }
    // TODO: better error message
    PyErr_SetString(PyExc_TypeError, "invalid item type");
    return -1;
}

static inline PyObject *VecNested_BoxItem(VecNested v, VecNestedBufItem item) {
    if (item.len < 0)
        Py_RETURN_NONE;
    VecNestedBufObject *v_buf = VEC_NESTED_BUF(v);
    VecNested_ItemXINCREF(v_buf, item);
    if (v_buf->depth > 1) {
        // Item is a nested vec
        VecNested iv = { .len = item.len,
            .items = (VecNestedBufItem *)item.items };
        return VecNested_Box(iv);
    } else {
        // Item is a non-nested vec
        size_t item_type = v_buf->item_type;
        if (item_type == VEC_ITEM_TYPE_I64) {
            VecI64 iv = { .len = item.len,
                .items = (int64_t *)item.items };
            return VecI64_Box(iv);
        } else if (item_type == VEC_ITEM_TYPE_U8) {
            VecU8 iv = { .len = item.len,
                .items = (uint8_t *)item.items };
            return VecU8_Box(iv);
        } else if (item_type == VEC_ITEM_TYPE_FLOAT) {
            VecFloat iv = { .len = item.len,
                .items = (double *)item.items };
            return VecFloat_Box(iv);
        } else if (item_type == VEC_ITEM_TYPE_I32) {
            VecI32 iv = { .len = item.len,
                .items = (int32_t *)item.items };
            return VecI32_Box(iv);
        } else if (item_type == VEC_ITEM_TYPE_I16) {
            VecI16 iv = { .len = item.len,
                .items = (int16_t *)item.items };
            return VecI16_Box(iv);
        } else if (item_type == VEC_ITEM_TYPE_BOOL) {
            VecBool iv = { .len = item.len,
                .items = (char *)item.items };
            return VecBool_Box(iv);
        } else {
            // Generic vec[t]
            VecT iv = { .len = item.len,
                .items = (PyObject **)item.items };
            return VecT_Box(iv, item_type);
        }
    }
}

// Growth helpers

static inline Py_ssize_t Vec_GrowCapacity(Py_ssize_t cap) {
    if (unlikely(cap > (PY_SSIZE_T_MAX - 1) / 2)) {
        // Allocation will fail at this size, but avoid overflow
        return PY_SSIZE_T_MAX;
    }
    return 2 * cap + 1;
}

static inline Py_ssize_t Vec_GrowCapacityTo(Py_ssize_t cap, Py_ssize_t min_cap) {
    while (cap < min_cap) {
        if (unlikely(cap > (PY_SSIZE_T_MAX - 1) / 2)) {
            cap = min_cap;
            break;
        }
        cap = 2 * cap + 1;
    }
    return cap;
}

// Misc helpers

PyObject *Vec_TypeToStr(size_t item_type, size_t depth);
PyObject *Vec_GenericRepr(PyObject *vec, size_t item_type, size_t depth, int verbose);
PyObject *Vec_GenericRichcompare(Py_ssize_t *len, PyObject **items,
                                  Py_ssize_t *other_len, PyObject **other_items,
                                  int op);
int Vec_GenericRemove(Py_ssize_t *len, PyObject **items, PyObject *item);
PyObject *Vec_GenericPopWrapper(Py_ssize_t *len, PyObject **items, PyObject *args);
PyObject *Vec_GenericPop(Py_ssize_t *len, PyObject **items, Py_ssize_t index);

#endif  // VEC_H_INCL
