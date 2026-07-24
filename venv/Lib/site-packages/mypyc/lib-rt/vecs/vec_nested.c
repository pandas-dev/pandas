// Implementation of nested vec[t], when t is a vec type.
//
// Examples of types supported:
//  - vec[vec[i64]]
//  - vec[vec[vec[str]]]
//  - vec[vec[str | None]]

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "librt_vecs.h"
#include "vecs_internal.h"

#define VEC_BUF(v) ((VecNestedBufObject *)((char *)(v).items - offsetof(VecNestedBufObject, items)))
#define VEC_CAP(v) (VEC_BUF(v)->ob_base.ob_size)
#define VEC_INCREF(v) do { if ((v).items) Py_INCREF(VEC_BUF(v)); } while (0)
#define VEC_DECREF(v) do { if ((v).items) Py_DECREF(VEC_BUF(v)); } while (0)

static inline VecNested vec_error() {
    VecNested v = { .len = -1 };
    return v;
}

static inline void vec_track_buffer(VecNested *vec) {
    PyObject_GC_Track(VEC_BUF(*vec));
}

static inline PyObject *box_vec_item_by_index(VecNested v, Py_ssize_t index) {
    return VecNested_BoxItem(v, v.items[index]);
}

// Alloc a partially initialized vec. If size > 0, caller *must* immediately initialize len,
// and items. Caller *must* also call vec_track_buffer on the returned vec but only
// after initializing the items.
static VecNested vec_alloc(Py_ssize_t size, size_t item_type, size_t depth) {
    VecNestedBufObject *buf = PyObject_GC_NewVar(VecNestedBufObject, &VecNestedBufType, size);
    if (buf == NULL)
        return vec_error();
    buf->item_type = item_type;
    buf->depth = depth;
    if (!Vec_IsMagicItemType(item_type))
        Py_INCREF(VEC_BUF_ITEM_TYPE(buf));
    VecNested res = { .items = buf->items };
    return res;
}

// Box a nested vec value, stealing 'vec'. On error, decref 'vec'.
PyObject *VecNested_Box(VecNested vec) {
    VecNestedObject *obj = PyObject_GC_New(VecNestedObject, &VecNestedType);
    if (obj == NULL) {
        VEC_DECREF(vec);
        return NULL;
    }
    obj->vec = vec;
    PyObject_GC_Track(obj);
    return (PyObject *)obj;
}

VecNested VecNested_Unbox(PyObject *obj, size_t item_type, size_t depth) {
    if (obj->ob_type == &VecNestedType) {
        VecNested result = ((VecNestedObject *)obj)->vec;
        VecNestedBufObject *buf = VEC_BUF(result);
        if (buf->item_type == item_type && buf->depth == depth) {
            VEC_INCREF(result);  // TODO: Should we borrow instead?
            return result;
        }
    }
    // TODO: Better error message, with name of type
    PyErr_SetString(PyExc_TypeError, "vec[t] expected");
    return vec_error();
}

VecNested VecNested_ConvertFromNested(VecNestedBufItem item) {
    return (VecNested) { item.len, (VecNestedBufItem *)item.items };
}

VecNested VecNested_New(Py_ssize_t size, Py_ssize_t cap, size_t item_type, size_t depth) {
    if (cap < 0) {
        PyErr_SetString(PyExc_ValueError, "capacity must not be negative");
        return vec_error();
    }
    if (cap < size)
        cap = size;
    VecNested vec = vec_alloc(cap, item_type, depth);
    if (VEC_IS_ERROR(vec))
        return vec;
    for (Py_ssize_t i = 0; i < cap; i++) {
        vec.items[i].len = -1;
        vec.items[i].items = 0;
    }
    vec.len = size;
    vec_track_buffer(&vec);
    return vec;
}

static PyObject *vec_repr(PyObject *self) {
    VecNested v = ((VecNestedObject *)self)->vec;
    VecNestedBufObject *buf = VEC_BUF(v);
    return Vec_GenericRepr(self, buf->item_type, buf->depth, 1);
}

static PyObject *vec_get_item(PyObject *o, Py_ssize_t i) {
    VecNested v = ((VecNestedObject *)o)->vec;
    if ((size_t)i < (size_t)v.len) {
        return box_vec_item_by_index(v, i);
    } else if ((size_t)i + (size_t)v.len < (size_t)v.len) {
        return box_vec_item_by_index(v, i + v.len);
    } else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
}

VecNested VecNested_Slice(VecNested vec, int64_t start, int64_t end) {
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
    VecNestedBufObject *vec_buf = VEC_BUF(vec);
    VecNested res = vec_alloc(slicelength, vec_buf->item_type, vec_buf->depth);
    if (VEC_IS_ERROR(res))
        return res;
    res.len = slicelength;
    for (Py_ssize_t i = 0; i < slicelength; i++) {
        VecNestedBufItem item = vec.items[start + i];
        VecNested_ItemXINCREF(vec_buf, item);
        res.items[i] = item;
    }
    vec_track_buffer(&res);
    return res;
}

static PyObject *vec_subscript(PyObject *self, PyObject *item) {
    VecNested vec = ((VecNestedObject *)self)->vec;
    if (PyIndex_Check(item)) {
        Py_ssize_t i = PyNumber_AsSsize_t(item, PyExc_IndexError);
        if (i == -1 && PyErr_Occurred())
            return NULL;
        if ((size_t)i < (size_t)vec.len) {
            return box_vec_item_by_index(vec, i);
        } else if ((size_t)i + (size_t)vec.len < (size_t)vec.len) {
            return box_vec_item_by_index(vec, i + vec.len);
        } else {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            return NULL;
        }
    } else if (PySlice_Check(item)) {
        Py_ssize_t start, stop, step;
        if (PySlice_Unpack(item, &start, &stop, &step) < 0)
            return NULL;
        Py_ssize_t slicelength = PySlice_AdjustIndices(vec.len, &start, &stop, step);
        VecNestedBufObject *vec_buf = VEC_BUF(vec);
        VecNested res = vec_alloc(slicelength, vec_buf->item_type, vec_buf->depth);
        if (VEC_IS_ERROR(res))
            return NULL;
        res.len = slicelength;
        Py_ssize_t j = start;
        for (Py_ssize_t i = 0; i < slicelength; i++) {
            VecNestedBufItem item = vec.items[j];
            VecNested_ItemXINCREF(vec_buf, item);
            res.items[i] = item;
            j += step;
        }
        vec_track_buffer(&res);
        return VecNested_Box(res);
    } else {
        PyErr_Format(PyExc_TypeError, "vec indices must be integers or slices, not %.100s",
                     item->ob_type->tp_name);
        return NULL;
    }
}

static int vec_ass_item(PyObject *self, Py_ssize_t i, PyObject *o) {
    VecNested v = ((VecNestedObject *)self)->vec;
    if ((size_t)i + (size_t)v.len < (size_t)v.len) {
        i += v.len;
    }
    if ((size_t)i < (size_t)v.len) {
        VecNestedBufItem item;
        if (VecNested_UnboxItem(v, o, &item) < 0)
            return -1;
        VecNestedBufObject *v_buf = VEC_BUF(v);
        VecNested_ItemXINCREF(v_buf, item);
        VecNested_ItemXDECREF(v_buf, v.items[i]);
        v.items[i] = item;
        return 0;
    } else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return -1;
    }
}

static int vec_contains(PyObject *self, PyObject *value) {
    VecNested v = ((VecNestedObject *)self)->vec;
    for (Py_ssize_t i = 0; i < v.len; i++) {
        PyObject *item = box_vec_item_by_index(v, i);
        if (item == NULL)
            return -1;
        if (item == value) {
            Py_DECREF(item);
            return 1;
        }
        int cmp = PyObject_RichCompareBool(item, value, Py_EQ);
        Py_DECREF(item);
        if (cmp != 0)
            return cmp;  // 1 if equal, -1 on error
    }
    return 0;
}

static PyObject *compare_vec_eq(VecNested x, VecNested y, int op) {
    int cmp = 1;
    PyObject *res;
    VecNestedBufObject *x_buf = VEC_BUF(x);
    VecNestedBufObject *y_buf = VEC_BUF(y);
    if (x.len != y.len
            || x_buf->item_type != y_buf->item_type
            || x_buf->depth != y_buf->depth) {
        cmp = 0;
    } else {
        for (Py_ssize_t i = 0; i < x.len; i++) {
            PyObject *x_item = box_vec_item_by_index(x, i);
            PyObject *y_item = box_vec_item_by_index(y, i);
            int itemcmp = PyObject_RichCompareBool(x_item, y_item, Py_EQ);
            Py_DECREF(x_item);
            Py_DECREF(y_item);
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
    Py_INCREF(res);
    return res;
}

PyObject *vec_richcompare(PyObject *self, PyObject *other, int op) {
    PyObject *res;
    if (op == Py_EQ || op == Py_NE) {
        if (other->ob_type != &VecNestedType) {
            res = op == Py_EQ ? Py_False : Py_True;
        } else {
            return compare_vec_eq(((VecNestedObject *)self)->vec,
                                        ((VecNestedObject *)other)->vec, op);
        }
    } else
        res = Py_NotImplemented;
    Py_INCREF(res);
    return res;
}

// Append item to 'vec', stealing 'vec'. Return 'vec' with item appended.
VecNested VecNested_Append(VecNested vec, VecNestedBufItem x) {
    Py_ssize_t cap = VEC_CAP(vec);
    VecNestedBufObject *vec_buf = VEC_BUF(vec);
    VecNested_ItemXINCREF(vec_buf, x);
    if (vec.len < cap) {
        // Slot may have duplicate ref from prior remove/pop
        VecNested_ItemXDECREF(vec_buf, vec.items[vec.len]);
        vec.items[vec.len] = x;
        vec.len++;
        return vec;
    } else {
        Py_ssize_t new_size = Vec_GrowCapacity(cap);
        // TODO: Avoid initializing to zero here
        VecNested new = vec_alloc(new_size, vec_buf->item_type, vec_buf->depth);
        if (VEC_IS_ERROR(new)) {
            VecNested_ItemXDECREF(vec_buf, x);
            // The input vec is being consumed/stolen by this function, so on error
            // we must decref it to avoid leaking the buffer.
            VEC_DECREF(vec);
            return new;
        }
        // Copy items to new vec.
        memcpy(new.items, vec.items, sizeof(VecNestedBufItem) * vec.len);
        // TODO: How to safely represent deleted items?
        memset(new.items + vec.len, 0, sizeof(VecNestedBufItem) * (new_size - vec.len));
        if (Py_REFCNT(vec_buf) > 1) {
            // Other references to old buffer exist; INCREF items in new buffer
            // so old buffer keeps valid references for aliases.
            for (Py_ssize_t i = 0; i < vec.len; i++)
                VecNested_ItemXINCREF(vec_buf, new.items[i]);
        } else {
            // No aliases; transfer ownership by clearing old buffer items.
            memset(vec.items, 0, sizeof(VecNestedBufItem) * vec.len);
        }
        new.items[vec.len] = x;
        new.len = vec.len + 1;
        vec_track_buffer(&new);
        VEC_DECREF(vec);
        return new;
    }
}

// Extend 'vec' with items from 'iterable', stealing 'vec'.
// Return extended 'vec', or error vec on failure.
VecNested VecNested_Extend(VecNested vec, PyObject *iterable) {
    if (VecNested_Check(iterable)) {
        VecNested src = ((VecNestedObject *)iterable)->vec;
        VecNestedBufObject *vec_buf = VEC_BUF(vec);
        VecNestedBufObject *src_buf = VEC_BUF(src);
        if (src_buf->item_type == vec_buf->item_type && src_buf->depth == vec_buf->depth) {
            return VecNested_ExtendVec(vec, src);
        }
    }

    PyObject *iter = PyObject_GetIter(iterable);
    if (iter == NULL) {
        VEC_DECREF(vec);
        return vec_error();
    }
    PyObject *item;
    while ((item = PyIter_Next(iter)) != NULL) {
        VecNestedBufItem vecitem;
        if (VecNested_UnboxItem(vec, item, &vecitem) < 0) {
            Py_DECREF(iter);
            VEC_DECREF(vec);
            Py_DECREF(item);
            return vec_error();
        }
        vec = VecNested_Append(vec, vecitem);
        Py_DECREF(item);
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
VecNested VecNested_ExtendVec(VecNested dst, VecNested src) {
    if (src.len == 0)
        return dst;
    if (src.len > PY_SSIZE_T_MAX - dst.len) {
        PyErr_NoMemory();
        VEC_DECREF(dst);
        return vec_error();
    }
    Py_ssize_t new_len = dst.len + src.len;
    // VecNested buf is never NULL (even for empty vecs), so no NULL guard needed
    Py_ssize_t cap = VEC_CAP(dst);
    VecNestedBufObject *dst_buf = VEC_BUF(dst);
    VecNestedBufObject *src_buf = VEC_BUF(src);
    if (new_len <= cap && dst.items != src.items) {
        // Fast path: enough capacity and no aliasing
        for (Py_ssize_t i = 0; i < src.len; i++) {
            VecNestedBufItem item = src.items[i];
            VecNested_ItemXINCREF(src_buf, item);
            // Slot may have duplicate ref from prior remove/pop
            VecNested_ItemXDECREF(dst_buf, dst.items[dst.len + i]);
            dst.items[dst.len + i] = item;
        }
        dst.len = new_len;
        return dst;
    }
    // Need to reallocate (or dst and src share a buffer)
    Py_ssize_t new_cap = Vec_GrowCapacityTo(cap, new_len);
    int aliased = dst.items == src.items;
    VecNested new = vec_alloc(new_cap, dst_buf->item_type, dst_buf->depth);
    if (VEC_IS_ERROR(new)) {
        VEC_DECREF(dst);
        return new;
    }
    if (aliased) {
        // dst and src share a buffer -- incref all items instead of
        // moving refs, to avoid mutating the shared buffer
        for (Py_ssize_t i = 0; i < dst.len; i++) {
            VecNested_ItemXINCREF(dst_buf, dst.items[i]);
            new.items[i] = dst.items[i];
        }
    } else {
        memcpy(new.items, dst.items, sizeof(VecNestedBufItem) * dst.len);
        if (Py_REFCNT(dst_buf) > 1) {
            for (Py_ssize_t i = 0; i < dst.len; i++)
                VecNested_ItemXINCREF(dst_buf, new.items[i]);
        } else {
            memset(dst.items, 0, sizeof(VecNestedBufItem) * dst.len);
        }
    }
    // Copy src items (incref each buf)
    for (Py_ssize_t i = 0; i < src.len; i++) {
        VecNestedBufItem item = src.items[i];
        VecNested_ItemXINCREF(src_buf, item);
        new.items[dst.len + i] = item;
    }
    memset(new.items + new_len, 0, sizeof(VecNestedBufItem) * (new_cap - new_len));
    new.len = new_len;
    vec_track_buffer(&new);
    VEC_DECREF(dst);
    return new;
}

// Remove item from 'vec', stealing 'vec'. Return 'vec' with item removed.
VecNested VecNested_Remove(VecNested self, VecNestedBufItem arg) {
    VecNestedBufItem *items = self.items;
    VecNestedBufObject *self_buf = VEC_BUF(self);

    PyObject *boxed_arg = VecNested_BoxItem(self, arg);
    if (boxed_arg == NULL) {
        // The input self is being consumed/stolen by this function, so on error
        // we must decref it to avoid leaking the buffer.
        VEC_DECREF(self);
        return vec_error();
    }

    for (Py_ssize_t i = 0; i < self.len; i++) {
        int match = 0;
        if (items[i].len == arg.len && items[i].items == arg.items)
            match = 1;
        else {
            PyObject *item = box_vec_item_by_index(self, i);
            if (item == NULL) {
                Py_DECREF(boxed_arg);
                // The input self is being consumed/stolen by this function, so on error
                // we must decref it to avoid leaking the buffer.
                VEC_DECREF(self);
                return vec_error();
            }
            int itemcmp = PyObject_RichCompareBool(item, boxed_arg, Py_EQ);
            Py_DECREF(item);
            if (itemcmp < 0) {
                Py_DECREF(boxed_arg);
                // The input self is being consumed/stolen by this function, so on error
                // we must decref it to avoid leaking the buffer.
                VEC_DECREF(self);
                return vec_error();
            }
            match = itemcmp;
        }
        if (match) {
            if (i < self.len - 1) {
                VecNested_ItemCLEAR(self_buf, &items[i]);
                for (; i < self.len - 1; i++) {
                    items[i] = items[i + 1];
                }
                VecNested_ItemXINCREF(self_buf, items[self.len - 1]);
            }
            self.len--;
            Py_DECREF(boxed_arg);
            // Return the stolen reference without INCREF
            return self;
        }
    }
    Py_DECREF(boxed_arg);
    PyErr_SetString(PyExc_ValueError, "vec.remove(x): x not in vec");
    // The input self is being consumed/stolen by this function, so on error
    // we must decref it to avoid leaking the buffer.
    VEC_DECREF(self);
    return vec_error();
}

// Pop item from 'vec', stealing 'vec'. Return struct with modified 'vec' and the popped item.
VecNestedPopResult VecNested_Pop(VecNested v, Py_ssize_t index) {
    VecNestedPopResult result;

    if (index < 0)
        index += v.len;

    if (index < 0 || index >= v.len) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        // The input v is being consumed/stolen by this function, so on error
        // we must decref it to avoid leaking the buffer.
        VEC_DECREF(v);
        result.f0 = vec_error();
        result.f1.len = 0;
        result.f1.items = NULL;
        return result;
    }

    VecNestedBufItem *items = v.items;
    VecNestedBufObject *v_buf = VEC_BUF(v);
    result.f1 = items[index];
    for (Py_ssize_t i = index; i < v.len - 1; i++)
        items[i] = items[i + 1];
    if (v.len > 0)
        VecNested_ItemXINCREF(v_buf, items[v.len - 1]);
    v.len--;
    // Return the stolen reference without INCREF
    result.f0 = v;
    return result;
}

static int
VecNested_traverse(VecNestedObject *self, visitproc visit, void *arg)
{
    if (self->vec.items)
        Py_VISIT(VEC_BUF(self->vec));
    return 0;
}

static int
VecNested_clear(VecNestedObject *self)
{
    if (self->vec.items) {
        Py_DECREF(VEC_BUF(self->vec));
        self->vec.items = NULL;
    }
    return 0;
}

static void
VecNested_dealloc(VecNestedObject *self)
{
    PyObject_GC_UnTrack(self);
    Py_TRASHCAN_BEGIN(self, VecNested_dealloc)
    VecNested_clear(self);
    Py_TYPE(self)->tp_free((PyObject *)self);
    Py_TRASHCAN_END
}

static int
VecNestedBuf_traverse(VecNestedBufObject *self, visitproc visit, void *arg)
{
    if (!Vec_IsMagicItemType(self->item_type))
        Py_VISIT(VEC_BUF_ITEM_TYPE(self));
    for (Py_ssize_t i = 0; i < VEC_BUF_SIZE(self); i++) {
        int ret = VecNested_ItemVISIT(self, self->items[i], visit, arg);
        if (ret)
            return ret;
    }
    return 0;
}

static inline int
VecNestedBuf_clear(VecNestedBufObject *self)
{
    for (Py_ssize_t i = 0; i < VEC_BUF_SIZE(self); i++) {
        VecNested_ItemCLEAR(self, &self->items[i]);
    }
    if (self->item_type && !Vec_IsMagicItemType(self->item_type)) {
        Py_DECREF(VEC_BUF_ITEM_TYPE(self));
        self->item_type = 0;
    }
    return 0;
}

static void
VecNestedBuf_dealloc(VecNestedBufObject *self)
{
    PyObject_GC_UnTrack(self);
    Py_TRASHCAN_BEGIN(self, VecNestedBuf_dealloc)
    VecNestedBuf_clear(self);
    Py_TYPE(self)->tp_free((PyObject *)self);
    Py_TRASHCAN_END
}

static Py_ssize_t vec_ext_length(PyObject *o) {
    return ((VecNestedObject *)o)->vec.len;
}

static PyMappingMethods VecNestedMapping = {
    .mp_length = vec_ext_length,
    .mp_subscript = vec_subscript,
};

static PySequenceMethods VecNestedSequence = {
    .sq_item = vec_get_item,
    .sq_ass_item = vec_ass_item,
    .sq_contains = vec_contains,
};

static PyMethodDef vec_methods[] = {
    {NULL, NULL, 0, NULL},  /* Sentinel */
};

// Iterator type for nested vecs

typedef struct {
    PyObject_HEAD
    VecNested vec;             // Unboxed vec (keeps buffer alive via items reference)
    Py_ssize_t index;          // Current iteration index
} VecNestedIterObject;

PyTypeObject VecNestedIterType;

static PyObject *VecNested_iter(PyObject *self) {
    VecNestedIterObject *it = PyObject_GC_New(VecNestedIterObject, &VecNestedIterType);
    if (it == NULL)
        return NULL;
    it->vec = ((VecNestedObject *)self)->vec;
    VEC_INCREF(it->vec);
    it->index = 0;
    PyObject_GC_Track(it);
    return (PyObject *)it;
}

static int
VecNestedIter_traverse(VecNestedIterObject *self, visitproc visit, void *arg)
{
    if (self->vec.items)
        Py_VISIT(VEC_BUF(self->vec));
    return 0;
}

static int
VecNestedIter_clear(VecNestedIterObject *self)
{
    if (self->vec.items) {
        Py_DECREF(VEC_BUF(self->vec));
        self->vec.items = NULL;
    }
    return 0;
}

static void VecNestedIter_dealloc(VecNestedIterObject *self) {
    PyObject_GC_UnTrack(self);
    VEC_DECREF(self->vec);
    PyObject_GC_Del(self);
}

static PyObject *VecNestedIter_next(VecNestedIterObject *self) {
    if (self->vec.items == NULL)
        return NULL;
    if (self->index < self->vec.len) {
        PyObject *item = box_vec_item_by_index(self->vec, self->index);
        if (item == NULL)
            return NULL;
        self->index++;
        return item;
    }
    VEC_DECREF(self->vec);
    self->vec.items = NULL;
    return NULL;  // StopIteration
}

static PyObject *VecNestedIter_len(VecNestedIterObject *self, PyObject *Py_UNUSED(ignored)) {
    if (self->vec.items == NULL)
        return PyLong_FromSsize_t(0);
    Py_ssize_t remaining = self->vec.len - self->index;
    if (remaining < 0)
        remaining = 0;
    return PyLong_FromSsize_t(remaining);
}

static PyMethodDef VecNestedIter_methods[] = {
    {"__length_hint__", (PyCFunction)VecNestedIter_len, METH_NOARGS, NULL},
    {NULL, NULL, 0, NULL},
};

PyTypeObject VecNestedIterType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "vec_nested_iterator",
    .tp_basicsize = sizeof(VecNestedIterObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
    .tp_traverse = (traverseproc)VecNestedIter_traverse,
    .tp_clear = (inquiry)VecNestedIter_clear,
    .tp_dealloc = (destructor)VecNestedIter_dealloc,
    .tp_iter = PyObject_SelfIter,
    .tp_iternext = (iternextfunc)VecNestedIter_next,
    .tp_methods = VecNestedIter_methods,
};

PyTypeObject VecNestedBufType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "vecbuf",
    .tp_doc = "Internal data buffer used by vec objects",
    .tp_basicsize = sizeof(VecNestedBufObject) - sizeof(VecNestedBufItem),
    .tp_itemsize = sizeof(VecNestedBufItem),
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
    .tp_traverse = (traverseproc)VecNestedBuf_traverse,
    //.tp_new = vecbuf_i64_new, //??
    .tp_free = PyObject_GC_Del,
    .tp_clear = (inquiry)VecNestedBuf_clear,
    .tp_dealloc = (destructor)VecNestedBuf_dealloc,
};

PyTypeObject VecNestedType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "vec",
    .tp_doc = "Mutable sequence-like container optimized for compilation with mypyc",
    .tp_basicsize = sizeof(VecNestedObject),
    .tp_itemsize = 0,
    .tp_base = &VecType,  // Inherit from base vec type for isinstance() support
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
    .tp_traverse = (traverseproc)VecNested_traverse,
    .tp_clear = (inquiry)VecNested_clear,
    .tp_dealloc = (destructor)VecNested_dealloc,
    //.tp_free = PyObject_GC_Del,
    .tp_repr = (reprfunc)vec_repr,
    .tp_iter = VecNested_iter,
    .tp_as_sequence = &VecNestedSequence,
    .tp_as_mapping = &VecNestedMapping,
    .tp_richcompare = vec_richcompare,
    .tp_methods = vec_methods,
    // TODO: free
};

PyObject *VecNested_FromIterable(size_t item_type, size_t depth, PyObject *iterable, int64_t cap) {
    VecNested v = vec_alloc(cap, item_type, depth);
    if (VEC_IS_ERROR(v))
        return NULL;
    if (cap > 0) {
        for (int64_t i = 0; i < cap; i++) {
            v.items[i].len = -1;
            v.items[i].items = 0;
        }
    }
    v.len = 0;
    vec_track_buffer(&v);

    PyObject *iter = PyObject_GetIter(iterable);
    if (iter == NULL) {
        VEC_DECREF(v);
        return NULL;
    }
    PyObject *item;
    while ((item = PyIter_Next(iter)) != NULL) {
        VecNestedBufItem vecitem;
        if (VecNested_UnboxItem(v, item, &vecitem) < 0) {
            Py_DECREF(iter);
            VEC_DECREF(v);
            Py_DECREF(item);
            return NULL;
        }
        v = VecNested_Append(v, vecitem);
        Py_DECREF(item);
        if (VEC_IS_ERROR(v)) {
            Py_DECREF(iter);
            VEC_DECREF(v);
            return NULL;
        }
    }
    Py_DECREF(iter);
    if (PyErr_Occurred()) {
        VEC_DECREF(v);
        return NULL;
    }
    return VecNested_Box(v);
}

VecNestedAPI Vec_NestedAPI = {
    &VecNestedType,
    &VecNestedBufType,
    VecNested_New,
    VecNested_Box,
    VecNested_Unbox,
    VecNested_ConvertFromNested,
    VecNested_Append,
    VecNested_Pop,
    VecNested_Remove,
    VecNested_Slice,
    VecNested_Extend,
    VecNested_ExtendVec,
};
