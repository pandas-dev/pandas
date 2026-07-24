// Definition of the mypyc.vecs module, which implements 'vec' and related functionality.
//
// NOTE: This is still experimental and work in progress.
//
// vec[t] is a sequence type that is optimized for use in mypyc-compiled code.
// Performance in interpreted code is a secondary concern. It's primarily optimized
// to provide fast item get/set operations and fast len(), especially for primitive
// value item types, such as 'i32', 'float' or 'bool'.
//
// Example:
//
//    from mypy_extensions import i64
//
//    v = vec[i64]()
//    v = append(v, 123)
//    print(v[0])
//    v[-1] = 234
//
// Key properties:
//
//  * In native code, a vec instance is an immutable value type with two fields
//    (length, item buffer).
//  * The item buffer can only store types of a specific type (enforced at runtime).
//  * There are specialized item buffer types for certain primitive value item types,
//    such as 'i32', 'float' and 'bool'. These use a packed value encoding.
//  * Since vec values are immutable (only the buffer is mutable), any operation
//    that changes the length (such as append) must return the updated vec value.
//  * The buffer often a has higher capacity than the length of a vec. Append
//    operations only allocate a new buffer if the capacity is exhausted, for O(1)
//    amortized time complexity for appends.
//  * Vec values can be boxed to PyObject *; this happens transparently and uses
//    a wrapper object.
//  * The item type can be a class type, or 't | None' if t is a class type. However,
//    optional value item types such as 'float | None' are not supported, and 'int'
//    is not valid as item type -- a native integer type such as `i64` or `u8` must
//    be used instead.
//  * Nested vecs are also supported (e.g. vec[vec[i32]]).
//  * All features are also available in interpreted code, but often at some
//    performance cost compared to using plain lists.
//
// Implementation summary:
//
//  * In a non-native context, instances are constructed via a temporary VecGenericAlias
//    object, which is the result of 'vec[<t>]' (indexing the 'vec' type).
//  * The module exports a C API via a capsule that mypyc compiled code uses for most
//    operations.
//  * There are multiple C structs that represent vec value objects with different item
//    types, such as VecI64 and VecT. The latter is for arbitrary reference objects
//    ('PyObject *' items internally).
//  * There are also multiple buffer types, such as VecI64Buf and VecTBuf.
//  * Similarly, there are multiple Python wrapper object types for boxed vecs, such as
//    VecI64Type and VecTType.
//  * An empty vec in compiled code can use a NULL pointer to represent the buffer,
//    saving an object allocation in a common use case.
//
// Motivation:
//  * A vec with value type items, such as vec[i64] or vec[vec[i32]], is very efficient
//    compared to a list object, since values are not boxed, there is no need for runtime
//    type checks on read, and index overflow checks are very fast due to the immutable
//    value encoding.
//  * A vec with a reference item type, such as vec[str], is often faster than a list,
//    since there is no need for a cast on item read and the index checks are faster
//    (among other things), but here the benefit is much smaller compared to using a
//    value item type.
//
// Summary of files:
//  * vec_t.c implements vec[t] and vec[t | None] for reference types
//  * vec_nested.c implements nested vecs, such as vec[vec[str]]
//  * vec_template.c is #included in multiples files (e.g. vec_i64.c) to produce
//    specialized code for vecs with value item types (e.g. vec[i64])
//  * This file defines the C extension and has some shared code

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "librt_vecs.h"

PyTypeObject *LibRTVecs_I64TypeObj;
PyTypeObject *LibRTVecs_I32TypeObj;
PyTypeObject *LibRTVecs_I16TypeObj;
PyTypeObject *LibRTVecs_U8TypeObj;


// vec generic alias
//
// Used for the result of vec[t] (indexing) in interpreted context that must preserve
// knowledge of 't'. These aren't real types. This only supports constructing instances.
typedef struct {
    PyObject_HEAD
    // Tagged pointer to PyTypeObject *, lowest bit set for optional item type
    // Can also be one of magic VEC_ITEM_TYPE_* constants
    size_t item_type;
    size_t depth;  // Number of nested VecNested or VecT types
} VecGenericAlias;

static PyObject *vec_generic_alias_call(PyObject *self, PyObject *args, PyObject *kw)
{
    static char *kwlist[] = {"", "capacity", NULL};
    PyObject *init = NULL;
    int64_t cap = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kw, "|OL:vec", kwlist, &init, &cap)) {
        return NULL;
    }
    if (cap < 0) {
        PyErr_SetString(PyExc_ValueError, "capacity must not be negative");
        return NULL;
    }
    VecGenericAlias *p = (VecGenericAlias *)self;
    if (p->depth == 0) {
        if (init == NULL) {
            VecT vec = VecT_New(0, cap, p->item_type);
            if (VEC_IS_ERROR(vec))
                return NULL;
            return VecT_Box(vec, p->item_type);
        } else {
            VecT vec = VecT_FromIterable(p->item_type, init, cap);
            if (VEC_IS_ERROR(vec))
                return NULL;
            return VecT_Box(vec, p->item_type);
        }
    } else {
        if (init == NULL) {
            VecNested vec = VecNested_New(0, cap, p->item_type, p->depth);
            if (VEC_IS_ERROR(vec))
                return NULL;
            return VecNested_Box(vec);
        } else {
            return VecNested_FromIterable(p->item_type, p->depth, init, cap);
        }
    }
}

static int
VecGenericAlias_traverse(VecGenericAlias *self, visitproc visit, void *arg)
{
    if (!Vec_IsMagicItemType(self->item_type))
        Py_VISIT((PyObject *)(self->item_type & ~1));
    return 0;
}

static void
VecGenericAlias_dealloc(VecGenericAlias *self)
{
    PyObject_GC_UnTrack(self);
    if (self->item_type && !Vec_IsMagicItemType(self->item_type)) {
        Py_DECREF((PyObject *)(self->item_type & ~1));
        self->item_type = 0;
    }
    PyObject_GC_Del(self);
}

PyObject *Vec_TypeToStr(size_t item_type, size_t depth) {
    PyObject *item = NULL;
    PyObject *result = NULL;

    if (depth == 0) {
        if ((item_type & ~1) == VEC_ITEM_TYPE_I64) {
            item = PyUnicode_FromFormat("i64");
        } else if ((item_type & ~1) == VEC_ITEM_TYPE_U8) {
            item = PyUnicode_FromFormat("u8");
        } else if ((item_type & ~1) == VEC_ITEM_TYPE_FLOAT) {
            item = PyUnicode_FromFormat("float");
        } else if ((item_type & ~1) == VEC_ITEM_TYPE_I32) {
            item = PyUnicode_FromFormat("i32");
        } else if ((item_type & ~1) == VEC_ITEM_TYPE_I16) {
            item = PyUnicode_FromFormat("i16");
        } else if ((item_type & ~1) == VEC_ITEM_TYPE_BOOL) {
            item = PyUnicode_FromFormat("bool");
        } else {
            item = PyObject_GetAttrString((PyObject *)(item_type & ~1), "__name__");
            if (item == NULL) {
                return NULL;
            }
            if (item_type & 1) {
                PyObject *optional_item = PyUnicode_FromFormat("%U | None", item);
                Py_DECREF(item);
                if (optional_item == NULL) {
                    return NULL;
                }
                item = optional_item;
            }
        }
    } else {
        item = Vec_TypeToStr(item_type, depth - 1);
    }

    if (item == NULL) {
        return NULL;
    }

    result = PyUnicode_FromFormat("vec[%U]", item);
    Py_DECREF(item);
    return result;
}

PyObject *VecGenericAlias_repr(PyObject *self) {
    PyObject *l = NULL;
    PyObject *prefix = NULL;
    PyObject *suffix = NULL;
    PyObject *sep = NULL;
    PyObject *type_str = NULL;
    PyObject *result = NULL;

    l = PyList_New(0);
    if (l == NULL) goto error;

    prefix = PyUnicode_FromString("<class_proxy '");
    if (prefix == NULL) goto error;

    suffix = PyUnicode_FromString("'>");
    if (suffix == NULL) goto error;

    sep = PyUnicode_FromString("");
    if (sep == NULL) goto error;

    if (PyList_Append(l, prefix) < 0) goto error;

    VecGenericAlias *v = (VecGenericAlias *)self;
    type_str = Vec_TypeToStr(v->item_type, v->depth);
    if (type_str == NULL) goto error;

    if (PyList_Append(l, type_str) < 0) goto error;
    if (PyList_Append(l, suffix) < 0) goto error;

    result = PyUnicode_Join(sep, l);

error:
    Py_XDECREF(l);
    Py_XDECREF(prefix);
    Py_XDECREF(suffix);
    Py_XDECREF(sep);
    Py_XDECREF(type_str);
    return result;
}

PyTypeObject VecGenericAliasType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "vec_generic_alias",
    .tp_basicsize = sizeof(VecGenericAlias),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
    .tp_call = vec_generic_alias_call,
    .tp_traverse = (traverseproc)VecGenericAlias_traverse,
    .tp_dealloc = (destructor)VecGenericAlias_dealloc,
    .tp_repr = (reprfunc)VecGenericAlias_repr,
};


// The 'vec' type
//
// This cannot be instantiated, and it's only used for isinstance and indexing: vec[T].
// All specialized vec types (VecI64Type, VecI32Type, etc.) inherit from this base type,
// so isinstance(v, vec) returns True for any specialized vec instance.

static PyObject *extract_optional_item(PyObject *item) {
    PyObject *args = PyObject_GetAttrString(item, "__args__");
    if (args == NULL) {
        PyErr_Clear();
        return NULL;
    }
    if (!PyTuple_CheckExact(args))
        goto error;
    if (PyTuple_GET_SIZE(args) != 2)
        goto error;
    PyObject *item0 = PyTuple_GET_ITEM(args, 0);
    PyObject *item1 = PyTuple_GET_ITEM(args, 1);
    if (item0 == (PyObject *)Py_None->ob_type) {
        Py_DECREF(args);
        return item1;
    } else if (item1 == (PyObject *)Py_None->ob_type) {
        Py_DECREF(args);
        return item0;
    }
  error:
    Py_DECREF(args);
    return NULL;
}

// Evaluation vec[<type>] in interpreted code and produce a generic alias object.
static PyObject *vec_class_getitem(PyObject *type, PyObject *item)
{
    if (item == (PyObject *)LibRTVecs_I64TypeObj) {
        Py_INCREF(&VecI64Type);
        return (PyObject *)&VecI64Type;
    } else if (item == (PyObject *)LibRTVecs_U8TypeObj) {
        Py_INCREF(&VecU8Type);
        return (PyObject *)&VecU8Type;
    } else if (item == (PyObject *)&PyFloat_Type) {
        Py_INCREF(&VecFloatType);
        return (PyObject *)&VecFloatType;
    } else if (item == (PyObject *)LibRTVecs_I32TypeObj) {
        Py_INCREF(&VecI32Type);
        return (PyObject *)&VecI32Type;
    } else if (item == (PyObject *)LibRTVecs_I16TypeObj) {
        Py_INCREF(&VecI16Type);
        return (PyObject *)&VecI16Type;
    } else if (item == (PyObject *)&PyBool_Type) {
        Py_INCREF(&VecBoolType);
        return (PyObject *)&VecBoolType;
    } else {
        size_t depth = 0;
        size_t item_type;
        int optional = 0;
        if (!PyObject_TypeCheck(item, &PyType_Type)) {
            PyObject *it = extract_optional_item(item);
            if (it != NULL) {
                optional = 1;
                item = it;
                // Check if this is a specialized value type that doesn't support None
                if (item == (PyObject *)LibRTVecs_I64TypeObj ||
                    item == (PyObject *)LibRTVecs_U8TypeObj ||
                    item == (PyObject *)&PyFloat_Type ||
                    item == (PyObject *)LibRTVecs_I32TypeObj ||
                    item == (PyObject *)LibRTVecs_I16TypeObj ||
                    item == (PyObject *)&PyBool_Type) {
                    PyErr_SetString(PyExc_TypeError,
                        "vec does not support optional value types (use vec[object] instead)");
                    return NULL;
                }
            }
            if (item->ob_type == &VecGenericAliasType) {
                if (optional) {
                    PyErr_SetString(PyExc_TypeError, "optional type not expected in vec[...]");
                    return NULL;
                }
                VecGenericAlias *p = (VecGenericAlias *)item;
                item_type = p->item_type;
                depth = p->depth + 1;
            } else if (!PyObject_TypeCheck(item, &PyType_Type)) {
                PyErr_SetString(PyExc_TypeError, "type object expected in vec[...]");
                return NULL;
            } else {
                item_type = (size_t)item | optional;
            }
        } else {
            if (item == (PyObject *)&VecI64Type) {
                item_type = VEC_ITEM_TYPE_I64;
                depth = 1;
            } else if (item == (PyObject *)&VecU8Type) {
                item_type = VEC_ITEM_TYPE_U8;
                depth = 1;
            } else if (item == (PyObject *)&VecFloatType) {
                item_type = VEC_ITEM_TYPE_FLOAT;
                depth = 1;
            } else if (item == (PyObject *)&VecI32Type) {
                item_type = VEC_ITEM_TYPE_I32;
                depth = 1;
            } else if (item == (PyObject *)&VecI16Type) {
                item_type = VEC_ITEM_TYPE_I16;
                depth = 1;
            } else if (item == (PyObject *)&VecBoolType) {
                item_type = VEC_ITEM_TYPE_BOOL;
                depth = 1;
            } else {
                item_type = (size_t)item;
            }
        }
        if (item == (PyObject *)&PyLong_Type) {
            PyErr_Format(PyExc_ValueError, "unsupported type in vec[%s] (use i64, i32, i16, or u8)",
                         ((PyTypeObject *)item)->tp_name);
            return NULL;
        }
        VecGenericAlias *p;
        p = PyObject_GC_New(VecGenericAlias, &VecGenericAliasType);
        if (p == NULL)
            return NULL;
        Py_INCREF(item);
        p->item_type = item_type;
        p->depth = depth;
        PyObject_GC_Track(p);
        return (PyObject *)p;
    }
}

static PyMethodDef vec_methods[] = {
    {"__class_getitem__", vec_class_getitem, METH_O|METH_CLASS, NULL},
    {NULL, NULL, 0, NULL},  /* Sentinel */
};

PyTypeObject VecType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "vec",
    .tp_basicsize = sizeof(VecBaseObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,  // No Py_TPFLAGS_BASETYPE - prevent Python subclassing
    .tp_methods = vec_methods,
};

// Common helpers

PyObject *Vec_GenericRepr(PyObject *vec, size_t item_type, size_t depth, int verbose) {
    PyObject *l = NULL;
    PyObject *prefix = NULL;
    PyObject *mid = NULL;
    PyObject *suffix = NULL;
    PyObject *sep = NULL;
    PyObject *comma = NULL;
    PyObject *result = NULL;

    l = PyList_New(0);
    if (l == NULL) goto error;

    sep = PyUnicode_FromString("");
    if (sep == NULL) goto error;

    comma = PyUnicode_FromString(", ");
    if (comma == NULL) goto error;

    if (verbose) {
        prefix = Vec_TypeToStr(item_type, depth);
        if (prefix == NULL) goto error;

        mid = PyUnicode_FromString("([");
        if (mid == NULL) goto error;

        suffix = PyUnicode_FromString("])");
        if (suffix == NULL) goto error;
    } else {
        prefix = PyUnicode_FromString("");
        if (prefix == NULL) goto error;

        mid = PyUnicode_FromString("[");
        if (mid == NULL) goto error;

        suffix = PyUnicode_FromString("]");
        if (suffix == NULL) goto error;
    }

    if (PyList_Append(l, prefix) < 0) goto error;
    if (PyList_Append(l, mid) < 0) goto error;

    Py_ssize_t len = PyObject_Length(vec);
    if (len < 0) goto error;

    for (Py_ssize_t i = 0; i < len; i++) {
        PyObject *it = PySequence_GetItem(vec, i);
        if (it == NULL) goto error;

        PyObject *r;
        if (depth == 0 || it == Py_None) {
            r = PyObject_Repr(it);
        } else {
            r = Vec_GenericRepr(it, item_type, depth - 1, 0);
        }
        Py_DECREF(it);

        if (r == NULL) goto error;

        if (PyList_Append(l, r) < 0) {
            Py_DECREF(r);
            goto error;
        }
        Py_DECREF(r);

        if (i + 1 < len) {
            if (PyList_Append(l, comma) < 0) goto error;
        }
    }

    if (PyList_Append(l, suffix) < 0) goto error;

    result = PyUnicode_Join(sep, l);

error:
    Py_XDECREF(l);
    Py_XDECREF(prefix);
    Py_XDECREF(mid);
    Py_XDECREF(suffix);
    Py_XDECREF(sep);
    Py_XDECREF(comma);
    return result;
}

// Generic comparison implementation for vecs with PyObject * items.
PyObject *Vec_GenericRichcompare(Py_ssize_t *len, PyObject **items,
                                  Py_ssize_t *other_len, PyObject **other_items,
                                  int op) {
    int cmp = 1;
    PyObject *res;
    if (op == Py_EQ || op == Py_NE) {
        if (*len != *other_len) {
            cmp = 0;
        } else {
            for (Py_ssize_t i = 0; i < *len && i < *other_len; i++) {
                int itemcmp = PyObject_RichCompareBool(items[i], other_items[i], Py_EQ);
                if (itemcmp < 0)
                    return NULL;
                if (!itemcmp) {
                    cmp = 0;
                    break;
                }
            }
        }
        if (op == Py_NE)
            cmp = cmp ^ 1;
        res = cmp ? Py_True : Py_False;
    } else
        res = Py_NotImplemented;
    Py_INCREF(res);
    return res;
}

int Vec_GenericRemove(Py_ssize_t *len, PyObject **items, PyObject *item) {
    for (Py_ssize_t i = 0; i < *len; i++) {
        int match = 0;
        if (items[i] == item)
            match = 1;
        else {
            int itemcmp = PyObject_RichCompareBool(items[i], item, Py_EQ);
            if (itemcmp < 0)
                return 0;
            match = itemcmp;
        }
        if (match) {
            Py_CLEAR(items[i]);
            for (; i < *len - 1; i++) {
                items[i] = items[i + 1];
            }
            (*len)--;
            return 1;
        }
    }
    PyErr_SetString(PyExc_ValueError, "vec.remove(x): x not in vec");
    return 0;
}

PyObject *Vec_GenericPop(Py_ssize_t *len, PyObject **items, Py_ssize_t index) {
    if (index < 0)
        index += *len;

    if (index < 0 || index >= *len) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }

    PyObject *item = items[index];
    for (Py_ssize_t i = index; i < *len - 1; i++)
        items[i] = items[i + 1];

    (*len)--;
    return item;
}

PyObject *Vec_GenericPopWrapper(Py_ssize_t *len, PyObject **items, PyObject *args) {
    Py_ssize_t index = -1;
    if (!PyArg_ParseTuple(args, "|n:pop", &index))
        return NULL;

    return Vec_GenericPop(len, items, index);
}

// Module-level functions

static PyObject *vec_append(PyObject *self, PyObject *args)
{
    PyObject *vec;
    PyObject *item;

    if (!PyArg_ParseTuple(args, "OO", &vec, &item))
        return NULL;

    if (VecI64_Check(vec)) {
        int64_t x = VecI64_UnboxItem(item);
        if (VecI64_IsUnboxError(x)) {
            return NULL;
        }
        VecI64 v = ((VecI64Object *)vec)->vec;
        VEC_I64_INCREF(v);
        v = VecI64_Append(v, x);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecI64_Box(v);
    } else if (VecU8_Check(vec)) {
        uint8_t x = VecU8_UnboxItem(item);
        if (VecU8_IsUnboxError(x)) {
            return NULL;
        }
        VecU8 v = ((VecU8Object *)vec)->vec;
        VEC_U8_INCREF(v);
        v = VecU8_Append(v, x);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecU8_Box(v);
    } else if (VecFloat_Check(vec)) {
        double x = VecFloat_UnboxItem(item);
        if (VecFloat_IsUnboxError(x)) {
            return NULL;
        }
        VecFloat v = ((VecFloatObject *)vec)->vec;
        VEC_FLOAT_INCREF(v);
        v = VecFloat_Append(v, x);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecFloat_Box(v);
    } else if (VecI32_Check(vec)) {
        int32_t x = VecI32_UnboxItem(item);
        if (VecI32_IsUnboxError(x)) {
            return NULL;
        }
        VecI32 v = ((VecI32Object *)vec)->vec;
        VEC_I32_INCREF(v);
        v = VecI32_Append(v, x);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecI32_Box(v);
    } else if (VecI16_Check(vec)) {
        int16_t x = VecI16_UnboxItem(item);
        if (VecI16_IsUnboxError(x)) {
            return NULL;
        }
        VecI16 v = ((VecI16Object *)vec)->vec;
        VEC_I16_INCREF(v);
        v = VecI16_Append(v, x);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecI16_Box(v);
    } else if (VecBool_Check(vec)) {
        double x = VecBool_UnboxItem(item);
        if (VecBool_IsUnboxError(x)) {
            return NULL;
        }
        VecBool v = ((VecBoolObject *)vec)->vec;
        VEC_BOOL_INCREF(v);
        v = VecBool_Append(v, x);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecBool_Box(v);
    } else if (VecT_Check(vec)) {
        VecT v = ((VecTObject *)vec)->vec;
        if (!VecT_ItemCheck(v, item, VEC_T_BUF(v)->item_type)) {
            return NULL;
        }
        VEC_T_INCREF(v);
        v = VecT_Append(v, item, VEC_T_BUF(v)->item_type);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecT_Box(v, VEC_T_BUF(v)->item_type);
    } else if (VecNested_Check(vec)) {
        VecNested v = ((VecNestedObject *)vec)->vec;
        VecNestedBufItem vecitem;
        if (VecNested_UnboxItem(v, item, &vecitem) < 0)
            return NULL;
        VEC_NESTED_INCREF(v);
        v = VecNested_Append(v, vecitem);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecNested_Box(v);
    } else {
        PyErr_Format(PyExc_TypeError, "vec argument expected, got %.100s",
                     Py_TYPE(vec)->tp_name);
        return NULL;
    }
}

static PyObject *vec_remove(PyObject *self, PyObject *args)
{
    PyObject *vec;
    PyObject *item;

    if (!PyArg_ParseTuple(args, "OO", &vec, &item))
        return NULL;

    if (VecI64_Check(vec)) {
        int64_t x = VecI64_UnboxItem(item);
        if (VecI64_IsUnboxError(x)) {
            return NULL;
        }
        VecI64 v = ((VecI64Object *)vec)->vec;
        VEC_I64_INCREF(v);
        v = VecI64_Remove(v, x);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecI64_Box(v);
    } else if (VecU8_Check(vec)) {
        uint8_t x = VecU8_UnboxItem(item);
        if (VecU8_IsUnboxError(x)) {
            return NULL;
        }
        VecU8 v = ((VecU8Object *)vec)->vec;
        VEC_U8_INCREF(v);
        v = VecU8_Remove(v, x);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecU8_Box(v);
    } else if (VecFloat_Check(vec)) {
        double x = VecFloat_UnboxItem(item);
        if (VecFloat_IsUnboxError(x)) {
            return NULL;
        }
        VecFloat v = ((VecFloatObject *)vec)->vec;
        VEC_FLOAT_INCREF(v);
        v = VecFloat_Remove(v, x);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecFloat_Box(v);
    } else if (VecI32_Check(vec)) {
        int32_t x = VecI32_UnboxItem(item);
        if (VecI32_IsUnboxError(x)) {
            return NULL;
        }
        VecI32 v = ((VecI32Object *)vec)->vec;
        VEC_I32_INCREF(v);
        v = VecI32_Remove(v, x);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecI32_Box(v);
    } else if (VecI16_Check(vec)) {
        int16_t x = VecI16_UnboxItem(item);
        if (VecI16_IsUnboxError(x)) {
            return NULL;
        }
        VecI16 v = ((VecI16Object *)vec)->vec;
        VEC_I16_INCREF(v);
        v = VecI16_Remove(v, x);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecI16_Box(v);
    } else if (VecBool_Check(vec)) {
        char x = VecBool_UnboxItem(item);
        if (VecBool_IsUnboxError(x)) {
            return NULL;
        }
        VecBool v = ((VecBoolObject *)vec)->vec;
        VEC_BOOL_INCREF(v);
        v = VecBool_Remove(v, x);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecBool_Box(v);
    } else if (VecT_Check(vec)) {
        VecT v = ((VecTObject *)vec)->vec;
        if (!VecT_ItemCheck(v, item, VEC_T_BUF(v)->item_type)) {
            return NULL;
        }
        VEC_T_INCREF(v);
        v = VecT_Remove(v, item);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecT_Box(v, VEC_T_BUF(v)->item_type);
    } else if (VecNested_Check(vec)) {
        VecNested v = ((VecNestedObject *)vec)->vec;
        VecNestedBufItem vecitem;
        if (VecNested_UnboxItem(v, item, &vecitem) < 0)
            return NULL;
        VEC_NESTED_INCREF(v);
        v = VecNested_Remove(v, vecitem);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecNested_Box(v);
    } else {
        PyErr_Format(PyExc_TypeError, "vec argument expected, got %.100s",
                     Py_TYPE(vec)->tp_name);
        return NULL;
    }
}

static PyObject *vec_pop(PyObject *self, PyObject *args)
{
    PyObject *vec;
    Py_ssize_t index = -1;

    if (!PyArg_ParseTuple(args, "O|n:pop", &vec, &index))
        return NULL;

    PyObject *result_item0;
    PyObject *result_item1;

    if (VecI64_Check(vec)) {
        VecI64 v = ((VecI64Object *)vec)->vec;
        VEC_I64_INCREF(v);
        VecI64PopResult r;
        r = VecI64_Pop(v, index);
        if (VEC_IS_ERROR(r.f0))
            return NULL;
        if ((result_item0 = VecI64_Box(r.f0)) == NULL)
            return NULL;
        result_item1 = VecI64_BoxItem(r.f1);
    } else if (VecU8_Check(vec)) {
        VecU8 v = ((VecU8Object *)vec)->vec;
        VEC_U8_INCREF(v);
        VecU8PopResult r;
        r = VecU8_Pop(v, index);
        if (VEC_IS_ERROR(r.f0))
            return NULL;
        if ((result_item0 = VecU8_Box(r.f0)) == NULL)
            return NULL;
        result_item1 = VecU8_BoxItem(r.f1);
    } else if (VecFloat_Check(vec)) {
        VecFloat v = ((VecFloatObject *)vec)->vec;
        VEC_FLOAT_INCREF(v);
        VecFloatPopResult r;
        r = VecFloat_Pop(v, index);
        if (VEC_IS_ERROR(r.f0))
            return NULL;
        if ((result_item0 = VecFloat_Box(r.f0)) == NULL)
            return NULL;
        result_item1 = VecFloat_BoxItem(r.f1);
    } else if (VecI32_Check(vec)) {
        VecI32 v = ((VecI32Object *)vec)->vec;
        VEC_I32_INCREF(v);
        VecI32PopResult r;
        r = VecI32_Pop(v, index);
        if (VEC_IS_ERROR(r.f0))
            return NULL;
        if ((result_item0 = VecI32_Box(r.f0)) == NULL)
            return NULL;
        result_item1 = VecI32_BoxItem(r.f1);
    } else if (VecI16_Check(vec)) {
        VecI16 v = ((VecI16Object *)vec)->vec;
        VEC_I16_INCREF(v);
        VecI16PopResult r;
        r = VecI16_Pop(v, index);
        if (VEC_IS_ERROR(r.f0))
            return NULL;
        if ((result_item0 = VecI16_Box(r.f0)) == NULL)
            return NULL;
        result_item1 = VecI16_BoxItem(r.f1);
    } else if (VecBool_Check(vec)) {
        VecBool v = ((VecBoolObject *)vec)->vec;
        VEC_BOOL_INCREF(v);
        VecBoolPopResult r;
        r = VecBool_Pop(v, index);
        if (VEC_IS_ERROR(r.f0))
            return NULL;
        if ((result_item0 = VecBool_Box(r.f0)) == NULL)
            return NULL;
        result_item1 = VecBool_BoxItem(r.f1);
    } else if (VecT_Check(vec)) {
        VecT v = ((VecTObject *)vec)->vec;
        VEC_T_INCREF(v);
        VecTPopResult r;
        r = VecT_Pop(v, index);
        if (VEC_IS_ERROR(r.f0))
            return NULL;
        result_item0 = VecT_Box(r.f0, VEC_T_BUF(v)->item_type);
        if (result_item0 == NULL) {
            Py_DECREF(r.f1);
            return NULL;
        }
        result_item1 = r.f1;
    } else if (VecNested_Check(vec)) {
        VecNested v = ((VecNestedObject *)vec)->vec;
        VEC_NESTED_INCREF(v);
        VecNestedPopResult r;
        r = VecNested_Pop(v, index);
        if (VEC_IS_ERROR(r.f0))
            return NULL;
        PyObject *popped_buf = VecNested_ItemBuf(VEC_NESTED_BUF(r.f0), r.f1);
        result_item0 = VecNested_Box(r.f0);
        if (result_item0 == NULL) {
            Py_XDECREF(popped_buf);
            return NULL;
        }
        result_item1 = VecNested_BoxItem(r.f0, r.f1);
        if (result_item1 == NULL) {
            Py_DECREF(result_item0);
            Py_XDECREF(popped_buf);
            return NULL;
        }
        Py_XDECREF(popped_buf);
    } else {
        PyErr_Format(PyExc_TypeError, "vec argument expected, got %.100s",
                     Py_TYPE(vec)->tp_name);
        return NULL;
    }

    if (result_item1 == NULL) {
        Py_DECREF(result_item0);
        return NULL;
    }

    PyObject *res = PyTuple_New(2);
    if (res == NULL) {
        Py_DECREF(result_item0);
        Py_DECREF(result_item1);
        return NULL;
    }

    PyTuple_SET_ITEM(res, 0, result_item0);
    PyTuple_SET_ITEM(res, 1, result_item1);
    return res;
}

static PyObject *vec_extend(PyObject *self, PyObject *args)
{
    PyObject *vec;
    PyObject *iterable;

    if (!PyArg_ParseTuple(args, "OO:extend", &vec, &iterable))
        return NULL;

    if (VecI64_Check(vec)) {
        VecI64 v = ((VecI64Object *)vec)->vec;
        VEC_I64_INCREF(v);
        v = VecI64_Extend(v, iterable);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecI64_Box(v);
    } else if (VecU8_Check(vec)) {
        VecU8 v = ((VecU8Object *)vec)->vec;
        VEC_U8_INCREF(v);
        v = VecU8_Extend(v, iterable);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecU8_Box(v);
    } else if (VecFloat_Check(vec)) {
        VecFloat v = ((VecFloatObject *)vec)->vec;
        VEC_FLOAT_INCREF(v);
        v = VecFloat_Extend(v, iterable);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecFloat_Box(v);
    } else if (VecI32_Check(vec)) {
        VecI32 v = ((VecI32Object *)vec)->vec;
        VEC_I32_INCREF(v);
        v = VecI32_Extend(v, iterable);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecI32_Box(v);
    } else if (VecI16_Check(vec)) {
        VecI16 v = ((VecI16Object *)vec)->vec;
        VEC_I16_INCREF(v);
        v = VecI16_Extend(v, iterable);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecI16_Box(v);
    } else if (VecBool_Check(vec)) {
        VecBool v = ((VecBoolObject *)vec)->vec;
        VEC_BOOL_INCREF(v);
        v = VecBool_Extend(v, iterable);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecBool_Box(v);
    } else if (VecT_Check(vec)) {
        VecT v = ((VecTObject *)vec)->vec;
        size_t item_type = VEC_T_BUF(v)->item_type;
        VEC_T_INCREF(v);
        v = VecT_Extend(v, iterable, item_type);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecT_Box(v, item_type);
    } else if (VecNested_Check(vec)) {
        VecNested v = ((VecNestedObject *)vec)->vec;
        VEC_NESTED_INCREF(v);
        v = VecNested_Extend(v, iterable);
        if (VEC_IS_ERROR(v))
            return NULL;
        return VecNested_Box(v);
    } else {
        PyErr_Format(PyExc_TypeError, "vec argument expected, got %.100s",
                     Py_TYPE(vec)->tp_name);
        return NULL;
    }
}

// Return the base VecType for isinstance checks
static PyTypeObject *get_vec_type(void) {
    return &VecType;
}

static int
vecs_abi_version(void) {
    return LIBRT_VECS_ABI_VERSION;
}

static int
vecs_api_version(void) {
    return LIBRT_VECS_API_VERSION;
}

static VecCapsule Capsule = {
    vecs_abi_version,
    vecs_api_version,
    &Vec_TAPI,
    &Vec_NestedAPI,
    &Vec_I64API,
    &Vec_I32API,
    &Vec_I16API,
    &Vec_U8API,
    &Vec_FloatAPI,
    &Vec_BoolAPI,
    get_vec_type,
};

static PyMethodDef VecsMethods[] = {
    {"append",  vec_append, METH_VARARGS, "Append a value to the end of a vec"},
    {"remove",  vec_remove, METH_VARARGS, "Remove first occurrence of value from a vec"},
    {"pop",  vec_pop, METH_VARARGS, "Remove and return vec item at index (default last)"},
    {"extend",  vec_extend, METH_VARARGS, "Extend a vec with items from an iterable"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static int
librt_vecs_module_exec(PyObject *m)
{
    PyObject *ext = PyImport_ImportModule("mypy_extensions");
    if (ext == NULL) {
        return -1;
    }

    LibRTVecs_I64TypeObj = (PyTypeObject *)PyObject_GetAttrString(ext, "i64");
    if (LibRTVecs_I64TypeObj == NULL) {
        return -1;
    }
    if (!PyType_Check(LibRTVecs_I64TypeObj)) {
        PyErr_SetString(PyExc_TypeError, "mypy_extensions.i64 is not a type");
        return -1;
    }
    LibRTVecs_I32TypeObj = (PyTypeObject *)PyObject_GetAttrString(ext, "i32");
    if (LibRTVecs_I32TypeObj == NULL) {
        return -1;
    }
    if (!PyType_Check(LibRTVecs_I32TypeObj)) {
        PyErr_SetString(PyExc_TypeError, "mypy_extensions.i32 is not a type");
        return -1;
    }
    LibRTVecs_I16TypeObj = (PyTypeObject *)PyObject_GetAttrString(ext, "i16");
    if (LibRTVecs_I16TypeObj == NULL) {
        return -1;
    }
    if (!PyType_Check(LibRTVecs_I16TypeObj)) {
        PyErr_SetString(PyExc_TypeError, "mypy_extensions.i16 is not a type");
        return -1;
    }
    LibRTVecs_U8TypeObj = (PyTypeObject *)PyObject_GetAttrString(ext, "u8");
    if (LibRTVecs_U8TypeObj == NULL) {
        return -1;
    }
    if (!PyType_Check(LibRTVecs_U8TypeObj)) {
        PyErr_SetString(PyExc_TypeError, "mypy_extensions.u8 is not a type");
        return -1;
    }

    if (PyType_Ready(&VecType) < 0)
        return -1;
    if (PyType_Ready(&VecGenericAliasType) < 0)
        return -1;

    if (PyType_Ready(&VecTType) < 0)
        return -1;
    if (PyType_Ready(&VecTBufType) < 0)
        return -1;
    if (PyType_Ready(&VecTIterType) < 0)
        return -1;

    if (PyType_Ready(&VecNestedType) < 0)
        return -1;
    if (PyType_Ready(&VecNestedBufType) < 0)
        return -1;
    if (PyType_Ready(&VecNestedIterType) < 0)
        return -1;

    if (PyType_Ready(&VecI64Type) < 0)
        return -1;
    if (PyType_Ready(&VecI64BufType) < 0)
        return -1;
    if (PyType_Ready(&VecI64IterType) < 0)
        return -1;
    if (PyType_Ready(&VecI32Type) < 0)
        return -1;
    if (PyType_Ready(&VecI32BufType) < 0)
        return -1;
    if (PyType_Ready(&VecI32IterType) < 0)
        return -1;
    if (PyType_Ready(&VecI16Type) < 0)
        return -1;
    if (PyType_Ready(&VecI16BufType) < 0)
        return -1;
    if (PyType_Ready(&VecI16IterType) < 0)
        return -1;
    if (PyType_Ready(&VecU8Type) < 0)
        return -1;
    if (PyType_Ready(&VecU8BufType) < 0)
        return -1;
    if (PyType_Ready(&VecU8IterType) < 0)
        return -1;
    if (PyType_Ready(&VecFloatType) < 0)
        return -1;
    if (PyType_Ready(&VecFloatBufType) < 0)
        return -1;
    if (PyType_Ready(&VecFloatIterType) < 0)
        return -1;
    if (PyType_Ready(&VecBoolType) < 0)
        return -1;
    if (PyType_Ready(&VecBoolBufType) < 0)
        return -1;
    if (PyType_Ready(&VecBoolIterType) < 0)
        return -1;

    Py_INCREF(&VecType);
    if (PyModule_AddObject(m, "vec", (PyObject *)&VecType) < 0) {
        Py_DECREF(&VecType);
        return -1;
    }

    PyObject *c_api = PyCapsule_New(&Capsule, "librt.vecs._C_API", NULL);
    if (c_api == NULL)
        return -1;

    if (PyModule_AddObject(m, "_C_API", c_api) < 0) {
        Py_XDECREF(c_api);
        Py_DECREF(&VecType);
        return -1;
    }

    Py_DECREF(ext);
    return 0;
}

static PyModuleDef_Slot librt_vecs_module_slots[] = {
    {Py_mod_exec, librt_vecs_module_exec},
#ifdef Py_MOD_GIL_NOT_USED
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};

static PyModuleDef vecsmodule = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "vecs",
    .m_doc = "This module implements the vec type and related functionality.",
    .m_size = 0,
    .m_methods = VecsMethods,
    .m_slots = librt_vecs_module_slots,
};

PyMODINIT_FUNC
PyInit_vecs(void)
{
    return PyModuleDef_Init(&vecsmodule);
}
