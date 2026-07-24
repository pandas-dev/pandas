// NOTE: This file can't be compiled on its own, it must be #included
//       with certain #defines set, as described below.
//
// Implementation of a vec class specialized to a specific item type, such
// as vec[i64] or vec[float]. Assume that certain #defines are provided that
// provide all the item type specific definitions:
//
//   VEC             vec C type (e.g. VecI32)
//   VEC_TYPE        boxed vec type object (e.g. VecI32Type)
//   VEC_OBJECT      boxed Python object struct (e.g. VecI32Object)
//   BUF_OBJECT      buffer Python object struct (e.g. VecI32BufObject)
//   BUF_TYPE        buffer type object (e.g. VecI32BufType)
//   NAME(suffix)    macro to create prefixed name with given suffix (e.g. VecI32##suffix)
//   FUNC(suffix)    macro to create prefixed function name with suffix (e.g. VecI32_##suffix)
//   ITEM_TYPE_STR   vec item type as C string literal (e.g. "i32")
//   ITEM_TYPE_MAGIC integer constant corresponding to the item type (e.g. VEC_ITEM_TYPE_I32)
//   ITEM_C_TYPE     C type used for items (e.g. int32_t)
//   FEATURES        capsule API struct name (e.g. Vec_I32API)
//   BOX_ITEM        C function to box item (e.g. VecI32_BoxItem)
//   UNBOX_ITEM      C function to unbox item (e.g. VecI32_UnboxItem)
//   IS_UNBOX_ERROR  C function to check for unbox error (e.g. VecI32_IsUnboxError)

#ifndef VEC
#error "VEC must be defined"
#endif

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "librt_vecs.h"
#include "vecs_internal.h"
#include "mypyc_util.h"

#define VEC_BUF(v) ((BUF_OBJECT *)((char *)(v).items - offsetof(BUF_OBJECT, items)))
#define VEC_CAP(v) (VEC_BUF(v)->ob_base.ob_size)
#define VEC_INCREF(v) do { if ((v).items) Py_INCREF(VEC_BUF(v)); } while (0)
#define VEC_DECREF(v) do { if ((v).items) Py_DECREF(VEC_BUF(v)); } while (0)

inline static VEC vec_error() {
    VEC v = { .len = -1 };
    return v;
}

// Alloc a partially initialized vec. Caller *must* initialize len.
static VEC vec_alloc(Py_ssize_t size)
{
    BUF_OBJECT *buf;
    /* TODO: Check for overflow */
    if (size == 0) {
        buf = NULL;
    } else {
        buf = PyObject_NewVar(BUF_OBJECT, &BUF_TYPE, size);
        if (buf == NULL)
            return vec_error();
    }
    VEC res = { .items = (buf != NULL) ? buf->items : NULL };
    return res;
}

static void vec_dealloc(VEC_OBJECT *self) {
    if (self->vec.items) {
        Py_DECREF(VEC_BUF(self->vec));
        self->vec.items = NULL;
    }
    PyObject_Del(self);
}

// Box a vec[<itemtype>] value, stealing 'vec'. On error, decref 'vec'.
PyObject *FUNC(Box)(VEC vec) {
    VEC_OBJECT *obj = PyObject_New(VEC_OBJECT, &VEC_TYPE);
    if (obj == NULL) {
        VEC_DECREF(vec);
        return NULL;
    }
    obj->vec = vec;
    return (PyObject *)obj;
}

VEC FUNC(Unbox)(PyObject *obj) {
    if (obj->ob_type == &VEC_TYPE) {
        VEC result = ((VEC_OBJECT *)obj)->vec;
        VEC_INCREF(result);  // TODO: Should we borrow instead?
        return result;
    } else {
        PyErr_SetString(PyExc_TypeError, "vec[" ITEM_TYPE_STR "] expected");
        return vec_error();
    }
}

VEC FUNC(ConvertFromNested)(VecNestedBufItem item) {
    return (VEC) { item.len, (ITEM_C_TYPE *)item.items };
}

VEC FUNC(New)(Py_ssize_t size, Py_ssize_t cap) {
    if (cap < 0) {
        PyErr_SetString(PyExc_ValueError, "capacity must not be negative");
        return vec_error();
    }
    if (cap < size)
        cap = size;
    VEC vec = vec_alloc(cap);
    if (VEC_IS_ERROR(vec))
        return vec;
    for (Py_ssize_t i = 0; i < cap; i++) {
        vec.items[i] = 0;
    }
    vec.len = size;
    return vec;
}

#ifdef BUFFER_FORMAT_CHAR_OK
inline static int buffer_format_matches(const char *fmt) {
    char c = *fmt;
    if (c == '@' || c == '=') {
        c = fmt[1];
    }
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    else if (c == '<') { c = fmt[1]; }
    else if (c == '>' || c == '!') { return 0; }
#else
    else if (c == '>') { c = fmt[1]; }
    else if (c == '<' || c == '!') { return 0; }
#endif
    return c != '\0' && BUFFER_FORMAT_CHAR_OK(c);
}

// Try to get a compatible buffer view from 'obj'. Return 1 if successful
// (view is filled and caller must call PyBuffer_Release), 0 if the object
// doesn't support buffer protocol or the format doesn't match (no cleanup
// needed), or -1 on error.
inline static int vec_get_buffer(PyObject *obj, Py_buffer *view) {
    if (PyObject_GetBuffer(obj, view, PyBUF_C_CONTIGUOUS | PyBUF_FORMAT) != 0) {
        PyErr_Clear();
        return 0;
    }
    if (view->ndim == 1
        && view->itemsize == sizeof(ITEM_C_TYPE)
        && buffer_format_matches(view->format)) {
        return 1;
    }
    PyBuffer_Release(view);
    return 0;
}
#endif

static inline VEC vec_from_sequence(PyObject *seq, int64_t cap, const int is_list) {
    Py_ssize_t n = is_list ? PyList_GET_SIZE(seq) : PyTuple_GET_SIZE(seq);
    Py_ssize_t alloc_size = n > cap ? n : cap;
    VEC v = vec_alloc(alloc_size);
    if (VEC_IS_ERROR(v))
        return vec_error();
    for (Py_ssize_t i = 0; i < n; i++) {
        PyObject *item = is_list ? PyList_GET_ITEM(seq, i) : PyTuple_GET_ITEM(seq, i);
        ITEM_C_TYPE x = UNBOX_ITEM(item);
        if (IS_UNBOX_ERROR(x)) {
            VEC_DECREF(v);
            return vec_error();
        }
        v.items[i] = x;
    }
    v.len = n;
    return v;
}

VEC FUNC(FromIterable)(PyObject *iterable, int64_t cap) {
    if (cap < 0) {
        PyErr_SetString(PyExc_ValueError, "capacity must not be negative");
        return vec_error();
    }

    if (ITEM_TYPE_MAGIC == VEC_ITEM_TYPE_U8 && PyBytes_CheckExact(iterable)) {
        Py_ssize_t n = PyBytes_GET_SIZE(iterable);
        Py_ssize_t alloc_size = n > cap ? n : cap;
        VEC v = vec_alloc(alloc_size);
        if (VEC_IS_ERROR(v))
            return vec_error();
        if (n > 0)
            memcpy(v.items, PyBytes_AS_STRING(iterable), n);
        v.len = n;
        return v;
    }

#ifdef BUFFER_FORMAT_CHAR_OK
    Py_buffer view;
    int buf_ok = vec_get_buffer(iterable, &view);
    if (buf_ok < 0)
        return vec_error();
    if (buf_ok) {
        Py_ssize_t n = view.len / (Py_ssize_t)sizeof(ITEM_C_TYPE);
        Py_ssize_t alloc_size = n > cap ? n : cap;
        VEC v = vec_alloc(alloc_size);
        if (VEC_IS_ERROR(v)) {
            PyBuffer_Release(&view);
            return vec_error();
        }
        if (n > 0) {
            memcpy(v.items, view.buf, n * sizeof(ITEM_C_TYPE));
        }
        v.len = n;
        PyBuffer_Release(&view);
        return v;
    }
#endif

    if (PyList_CheckExact(iterable)) {
        return vec_from_sequence(iterable, cap, 1);
    } else if (PyTuple_CheckExact(iterable)) {
        return vec_from_sequence(iterable, cap, 0);
    }

    VEC v = vec_alloc(cap);
    if (VEC_IS_ERROR(v))
        return vec_error();
    if (cap > 0) {
        memset(v.items, 0, sizeof(ITEM_C_TYPE) * cap);
    }
    v.len = 0;

    PyObject *iter = PyObject_GetIter(iterable);
    if (iter == NULL) {
        VEC_DECREF(v);
        return vec_error();
    }
    PyObject *item;
    while ((item = PyIter_Next(iter)) != NULL) {
        ITEM_C_TYPE x = UNBOX_ITEM(item);
        Py_DECREF(item);
        if (IS_UNBOX_ERROR(x)) {
            Py_DECREF(iter);
            VEC_DECREF(v);
            return vec_error();
        }
        v = FUNC(Append)(v, x);
        if (VEC_IS_ERROR(v)) {
            Py_DECREF(iter);
            VEC_DECREF(v);
            return vec_error();
        }
    }
    Py_DECREF(iter);
    if (PyErr_Occurred()) {
        VEC_DECREF(v);
        return vec_error();
    }
    return v;
}

static PyObject *vec_new(PyTypeObject *self, PyObject *args, PyObject *kw) {
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
    if (init == NULL) {
        return FUNC(Box)(FUNC(New)(0, cap));
    } else {
        VEC v = FUNC(FromIterable)(init, cap);
        if (VEC_IS_ERROR(v))
            return NULL;
        return FUNC(Box)(v);
    }
}

static PyObject *vec_repr(PyObject *self) {
    return Vec_GenericRepr(self, ITEM_TYPE_MAGIC, 0, 1);
}

static PyObject *vec_get_item(PyObject *o, Py_ssize_t i) {
    VEC v = ((VEC_OBJECT *)o)->vec;
    if ((size_t)i < (size_t)v.len) {
        return BOX_ITEM(v.items[i]);
    } else if ((size_t)i + (size_t)v.len < (size_t)v.len) {
        return BOX_ITEM(v.items[i + v.len]);
    } else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
}

VEC FUNC(Slice)(VEC vec, int64_t start, int64_t end) {
    if (start < 0)
        start += vec.len;
    if (end < 0)
        end += vec.len;
    if (start < 0)
        start = 0;
    if (start >= vec.len)
        start = vec.len;
    if (end < start)
        end = start;
    if (end > vec.len)
        end = vec.len;
    int64_t slicelength = end - start;
    VEC res = vec_alloc(slicelength);
    if (VEC_IS_ERROR(res))
        return res;
    res.len = slicelength;
    for (Py_ssize_t i = 0; i < slicelength; i++)
        res.items[i] = vec.items[start + i];
    return res;
}

static PyObject *vec_subscript(PyObject *self, PyObject *item) {
    VEC vec = ((VEC_OBJECT *)self)->vec;
    if (PyIndex_Check(item)) {
        Py_ssize_t i = PyNumber_AsSsize_t(item, PyExc_IndexError);
        if (i == -1 && PyErr_Occurred())
            return NULL;
        if ((size_t)i < (size_t)vec.len) {
            return BOX_ITEM(vec.items[i]);
        } else if ((size_t)i + (size_t)vec.len < (size_t)vec.len) {
            return BOX_ITEM(vec.items[i + vec.len]);
        } else {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            return NULL;
        }
    } else if (PySlice_Check(item)) {
        Py_ssize_t start, stop, step;
        if (PySlice_Unpack(item, &start, &stop, &step) < 0)
            return NULL;
        Py_ssize_t slicelength = PySlice_AdjustIndices(vec.len, &start, &stop, step);
        VEC res = vec_alloc(slicelength);
        if (VEC_IS_ERROR(res))
            return NULL;
        res.len = slicelength;
        Py_ssize_t j = start;
        for (Py_ssize_t i = 0; i < slicelength; i++) {
            res.items[i] = vec.items[j];
            j += step;
        }
        return FUNC(Box)(res);
    } else {
        PyErr_Format(PyExc_TypeError, "vec indices must be integers or slices, not %.100s",
                     item->ob_type->tp_name);
        return NULL;
    }
}

static int vec_ass_item(PyObject *self, Py_ssize_t i, PyObject *o) {
    ITEM_C_TYPE x = UNBOX_ITEM(o);
    if (IS_UNBOX_ERROR(x))
        return -1;
    VEC v = ((VEC_OBJECT *)self)->vec;
    if ((size_t)i < (size_t)v.len) {
        v.items[i] = x;
        return 0;
    } else if ((size_t)i + (size_t)v.len < (size_t)v.len) {
        v.items[i + v.len] = x;
        return 0;
    } else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return -1;
    }
}

static int vec_contains(PyObject *self, PyObject *value) {
    ITEM_C_TYPE x = UNBOX_ITEM(value);
    if (unlikely(IS_UNBOX_ERROR(x))) {
        if (PyErr_Occurred())
            PyErr_Clear();
        // Fall back to boxed comparison (e.g. 2.0 == 2)
        VEC v = ((VEC_OBJECT *)self)->vec;
        for (Py_ssize_t i = 0; i < v.len; i++) {
            PyObject *boxed = BOX_ITEM(v.items[i]);
            if (boxed == NULL)
                return -1;
            int cmp = PyObject_RichCompareBool(boxed, value, Py_EQ);
            Py_DECREF(boxed);
            if (cmp != 0)
                return cmp;  // 1 if equal, -1 on error
        }
        return 0;
    }
    VEC v = ((VEC_OBJECT *)self)->vec;
    for (Py_ssize_t i = 0; i < v.len; i++) {
        if (v.items[i] == x)
            return 1;
    }
    return 0;
}

static Py_ssize_t vec_length(PyObject *o) {
    return ((VEC_OBJECT *)o)->vec.len;
}

static PyObject *vec_richcompare(PyObject *self, PyObject *other, int op) {
    int cmp = 1;
    PyObject *res;
    if (op == Py_EQ || op == Py_NE) {
        if (other->ob_type != &VEC_TYPE)
            cmp = 0;
        else {
            VEC x = ((VEC_OBJECT *)self)->vec;
            VEC y = ((VEC_OBJECT *)other)->vec;
            if (x.len != y.len) {
                cmp = 0;
            } else {
                for (Py_ssize_t i = 0; i < x.len; i++) {
                    if (x.items[i] != y.items[i]) {
                        cmp = 0;
                        break;
                    }
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

// Append item to 'vec', stealing 'vec'. Return 'vec' with item appended.
VEC FUNC(Append)(VEC vec, ITEM_C_TYPE x) {
    if (vec.items && vec.len < VEC_CAP(vec)) {
        vec.items[vec.len] = x;
        vec.len++;
        return vec;
    } else {
        Py_ssize_t cap = vec.items ? VEC_CAP(vec) : 0;
        Py_ssize_t new_size = Vec_GrowCapacity(cap);
        VEC new = vec_alloc(new_size);
        if (VEC_IS_ERROR(new)) {
            // The input v is being consumed/stolen by this function, so on error
            // we must decref it to avoid leaking the buffer.
            VEC_DECREF(vec);
            return vec_error();
        }
        new.len = vec.len + 1;
        if (vec.len > 0)
            memcpy(new.items, vec.items, sizeof(ITEM_C_TYPE) * vec.len);
        new.items[vec.len] = x;
        VEC_DECREF(vec);
        return new;
    }
}

inline static int vec_memory_overlaps(const void *p1, Py_ssize_t len1,
                                      const void *p2, Py_ssize_t len2) {
    if (len1 <= 0 || len2 <= 0)
        return 0;
    uintptr_t a = (uintptr_t)p1, b = (uintptr_t)p2;
    if (a <= b)
        return b - a < (uintptr_t)len1;
    return a - b < (uintptr_t)len2;
}

// Extend 'dst' by appending 'n' items from 'items', stealing 'dst'.
// Caller guarantees n > 0 and that 'items' remains valid for the call.
// If force_alloc is true, always allocate a new buffer even when dst has capacity.
inline static VEC vec_extend_items(
    VEC dst, const ITEM_C_TYPE *items, Py_ssize_t n, int force_alloc
) {
    if (unlikely(n > PY_SSIZE_T_MAX - dst.len)) {
        PyErr_NoMemory();
        VEC_DECREF(dst);
        return vec_error();
    }
    Py_ssize_t new_len = dst.len + n;
    Py_ssize_t cap = dst.items ? VEC_CAP(dst) : 0;
    if (!force_alloc && new_len <= cap) {
        memcpy(dst.items + dst.len, items, sizeof(ITEM_C_TYPE) * n);
        dst.len = new_len;
        return dst;
    }
    Py_ssize_t new_cap = Vec_GrowCapacityTo(cap, new_len);
    VEC new = vec_alloc(new_cap);
    if (VEC_IS_ERROR(new)) {
        VEC_DECREF(dst);
        return vec_error();
    }
    if (dst.len > 0)
        memcpy(new.items, dst.items, sizeof(ITEM_C_TYPE) * dst.len);
    memcpy(new.items + dst.len, items, sizeof(ITEM_C_TYPE) * n);
    new.len = new_len;
    VEC_DECREF(dst);
    return new;
}

// Extend 'vec' with items from 'iterable', stealing 'vec'.
// Return extended 'vec', or error vec on failure.
VEC FUNC(Extend)(VEC vec, PyObject *iterable) {
    if (Py_TYPE(iterable) == &VEC_TYPE) {
        return FUNC(ExtendVec)(vec, ((VEC_OBJECT *)iterable)->vec);
    }

    if (ITEM_TYPE_MAGIC == VEC_ITEM_TYPE_U8 && PyBytes_CheckExact(iterable)) {
        Py_ssize_t n = PyBytes_GET_SIZE(iterable);
        if (n > 0)
            return vec_extend_items(vec, (const ITEM_C_TYPE *)PyBytes_AS_STRING(iterable), n, 0);
        return vec;
    }

#ifdef BUFFER_FORMAT_CHAR_OK
    Py_buffer view;
    int buf_ok = vec_get_buffer(iterable, &view);
    if (buf_ok < 0) {
        VEC_DECREF(vec);
        return vec_error();
    }
    if (buf_ok) {
        Py_ssize_t n = view.len / (Py_ssize_t)sizeof(ITEM_C_TYPE);
        if (n > 0) {
            Py_ssize_t dst_bytes = n * (Py_ssize_t)sizeof(ITEM_C_TYPE);
            int force_alloc = vec.items != NULL
                && n <= VEC_CAP(vec) - vec.len
                && vec_memory_overlaps(view.buf, view.len,
                                       vec.items + vec.len, dst_bytes);
            vec = vec_extend_items(vec, (const ITEM_C_TYPE *)view.buf, n, force_alloc);
        }
        PyBuffer_Release(&view);
        return vec;
    }
#endif

    PyObject *iter = PyObject_GetIter(iterable);
    if (iter == NULL) {
        VEC_DECREF(vec);
        return vec_error();
    }
    PyObject *item;
    while ((item = PyIter_Next(iter)) != NULL) {
        ITEM_C_TYPE x = UNBOX_ITEM(item);
        Py_DECREF(item);
        if (IS_UNBOX_ERROR(x)) {
            Py_DECREF(iter);
            VEC_DECREF(vec);
            return vec_error();
        }
        vec = FUNC(Append)(vec, x);
        if (VEC_IS_ERROR(vec)) {
            Py_DECREF(iter);
            return vec_error();
        }
    }
    Py_DECREF(iter);
    if (PyErr_Occurred()) {
        VEC_DECREF(vec);
        return vec_error();
    }
    return vec;
}

// Extend 'dst' with items from 'src' vec, stealing 'dst', borrowing 'src'.
// Return extended vec, or error vec on failure.
VEC FUNC(ExtendVec)(VEC dst, VEC src) {
    if (src.len == 0)
        return dst;
    return vec_extend_items(dst, src.items, src.len, dst.items == src.items);
}

// Convert vec to list, stealing 'v'.
PyObject *FUNC(ToList)(VEC v) {
    Py_ssize_t n = v.len;
    PyObject *list = PyList_New(n);
    if (list == NULL) {
        VEC_DECREF(v);
        return NULL;
    }
    for (Py_ssize_t i = 0; i < n; i++) {
        PyObject *item = BOX_ITEM(v.items[i]);
        if (item == NULL) {
            Py_DECREF(list);
            VEC_DECREF(v);
            return NULL;
        }
        PyList_SET_ITEM(list, i, item);
    }
    VEC_DECREF(v);
    return list;
}

// Convert vec to tuple, stealing 'v'.
PyObject *FUNC(ToTuple)(VEC v) {
    Py_ssize_t n = v.len;
    PyObject *tuple = PyTuple_New(n);
    if (tuple == NULL) {
        VEC_DECREF(v);
        return NULL;
    }
    for (Py_ssize_t i = 0; i < n; i++) {
        PyObject *item = BOX_ITEM(v.items[i]);
        if (item == NULL) {
            Py_DECREF(tuple);
            VEC_DECREF(v);
            return NULL;
        }
        PyTuple_SET_ITEM(tuple, i, item);
    }
    VEC_DECREF(v);
    return tuple;
}

// Remove item from 'vec', stealing 'vec'. Return 'vec' with item removed.
VEC FUNC(Remove)(VEC v, ITEM_C_TYPE x) {
    for (Py_ssize_t i = 0; i < v.len; i++) {
        if (v.items[i] == x) {
            for (; i < v.len - 1; i++) {
                v.items[i] = v.items[i + 1];
            }
            v.len--;
            // Return the stolen reference without INCREF
            return v;
        }
    }
    PyErr_SetString(PyExc_ValueError, "vec.remove(x): x not in vec");
    // The input v is being consumed/stolen by this function, so on error
    // we must decref it to avoid leaking the buffer.
    VEC_DECREF(v);
    return vec_error();
}

// Pop item from 'vec', stealing 'vec'. Return struct with modified 'vec' and the popped item.
NAME(PopResult) FUNC(Pop)(VEC v, Py_ssize_t index) {
    NAME(PopResult) result;

    if (index < 0)
        index += v.len;

    if (index < 0 || index >= v.len) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        // The input v is being consumed/stolen by this function, so on error
        // we must decref it to avoid leaking the buffer.
        VEC_DECREF(v);
        result.f0 = vec_error();
        result.f1 = 0;
        return result;
    }

    result.f1 = v.items[index];
    for (Py_ssize_t i = index; i < v.len - 1; i++) {
        v.items[i] = v.items[i + 1];
    }

    v.len--;
    // Return the stolen reference without INCREF
    result.f0 = v;
    return result;
}

static PyMappingMethods vec_mapping_methods = {
    .mp_length = vec_length,
    .mp_subscript = vec_subscript,
};

#ifdef BUFFER_FORMAT
static int vec_getbuffer(VEC_OBJECT *self, Py_buffer *view, int flags) {
    if (view == NULL) {
        PyErr_SetString(PyExc_BufferError,
            "vec_getbuffer: view==NULL argument is obsolete");
        return -1;
    }
    if ((flags & PyBUF_WRITABLE) == PyBUF_WRITABLE) {
        PyErr_SetString(PyExc_BufferError, "Object is not writable");
        view->obj = NULL;
        return -1;
    }

    view->obj = (PyObject *)self;
    Py_INCREF(self);
    view->buf = (void *)self->vec.items;
    view->len = self->vec.len * (Py_ssize_t)sizeof(ITEM_C_TYPE);
    view->readonly = 1;
    view->itemsize = sizeof(ITEM_C_TYPE);
    view->format = NULL;
    if ((flags & PyBUF_FORMAT) == PyBUF_FORMAT)
        view->format = BUFFER_FORMAT;
    view->ndim = 1;
    view->shape = NULL;
    if ((flags & PyBUF_ND) == PyBUF_ND)
        view->shape = &self->vec.len;
    view->strides = NULL;
    if ((flags & PyBUF_STRIDES) == PyBUF_STRIDES)
        view->strides = &view->itemsize;
    view->suboffsets = NULL;
    view->internal = NULL;

    return 0;
}

static PyBufferProcs vec_buffer_procs = {
    .bf_getbuffer = (getbufferproc)vec_getbuffer,
};
#endif

static PySequenceMethods vec_sequence_methods = {
    .sq_item = vec_get_item,
    .sq_ass_item = vec_ass_item,
    .sq_contains = vec_contains,
};

static PyMethodDef vec_methods[] = {
    {NULL, NULL, 0, NULL},  /* Sentinel */
};

// Iterator type for specialized vec types

typedef struct {
    PyObject_HEAD
    VEC vec;              // Unboxed vec (keeps buffer alive via items reference)
    Py_ssize_t index;     // Current iteration index
} NAME(IterObject);

PyTypeObject NAME(IterType);

static PyObject *vec_iter(PyObject *self) {
    NAME(IterObject) *it = PyObject_New(NAME(IterObject), &NAME(IterType));
    if (it == NULL)
        return NULL;
    it->vec = ((VEC_OBJECT *)self)->vec;
    VEC_INCREF(it->vec);
    it->index = 0;
    return (PyObject *)it;
}

static void vec_iter_dealloc(NAME(IterObject) *self) {
    VEC_DECREF(self->vec);
    PyObject_Del(self);
}

static PyObject *vec_iter_next(NAME(IterObject) *self) {
    if (self->vec.items == NULL)
        return NULL;
    if (self->index < self->vec.len) {
        PyObject *item = BOX_ITEM(self->vec.items[self->index]);
        if (item == NULL)
            return NULL;
        self->index++;
        return item;
    }
    VEC_DECREF(self->vec);
    self->vec.items = NULL;
    return NULL;  // StopIteration
}

static PyObject *vec_iter_len(NAME(IterObject) *self, PyObject *Py_UNUSED(ignored)) {
    if (self->vec.items == NULL)
        return PyLong_FromSsize_t(0);
    Py_ssize_t remaining = self->vec.len - self->index;
    if (remaining < 0)
        remaining = 0;
    return PyLong_FromSsize_t(remaining);
}

static PyMethodDef vec_iter_methods[] = {
    {"__length_hint__", (PyCFunction)vec_iter_len, METH_NOARGS, NULL},
    {NULL, NULL, 0, NULL},
};

PyTypeObject NAME(IterType) = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "vec_iterator[" ITEM_TYPE_STR "]",
    .tp_basicsize = sizeof(NAME(IterObject)),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_dealloc = (destructor)vec_iter_dealloc,
    .tp_iter = PyObject_SelfIter,
    .tp_iternext = (iternextfunc)vec_iter_next,
    .tp_methods = vec_iter_methods,
};

PyTypeObject BUF_TYPE = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "vecbuf[" ITEM_TYPE_STR "]",
    .tp_doc = "Internal data buffer used by vec objects",
    .tp_basicsize = sizeof(BUF_OBJECT) - sizeof(ITEM_C_TYPE),
    .tp_itemsize = sizeof(ITEM_C_TYPE),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    //.tp_new = ??
    .tp_free = PyObject_Del,
};

PyTypeObject VEC_TYPE = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "vec[" ITEM_TYPE_STR "]",
    .tp_doc = "Mutable sequence-like container optimized for compilation with mypyc",
    .tp_basicsize = sizeof(VEC_OBJECT),
    .tp_itemsize = 0,
    .tp_base = &VecType,  // Inherit from base vec type for isinstance() support
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = vec_new,
    //.tp_free = PyObject_Del,
    .tp_dealloc = (destructor)vec_dealloc,
    .tp_repr = (reprfunc)vec_repr,
    .tp_iter = vec_iter,
    .tp_as_sequence = &vec_sequence_methods,
    .tp_as_mapping = &vec_mapping_methods,
#ifdef BUFFER_FORMAT
    .tp_as_buffer = &vec_buffer_procs,
#endif
    .tp_richcompare = vec_richcompare,
    .tp_methods = vec_methods,
};

NAME(API) FEATURES = {
    &VEC_TYPE,
    &BUF_TYPE,
    FUNC(New),
    FUNC(Box),
    FUNC(Unbox),
    FUNC(ConvertFromNested),
    FUNC(Append),
    FUNC(Pop),
    FUNC(Remove),
    FUNC(Slice),
    FUNC(FromIterable),
    FUNC(Extend),
    FUNC(ExtendVec),
    FUNC(ToList),
    FUNC(ToTuple),
};

#undef VEC_BUF
#undef VEC_CAP
#undef VEC_INCREF
#undef VEC_DECREF
