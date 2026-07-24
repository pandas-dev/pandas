#include "pythoncapi_compat.h"

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdint.h>
#include "CPy.h"
#include "librt_strings.h"

#define CPY_BOOL_ERROR 2

//
// BytesWriter
//

#define _WRITE_BYTES(data, type, v) \
    do { \
       *(type *)(((BytesWriterObject *)data)->buf + ((BytesWriterObject *)data)->len) = v; \
       ((BytesWriterObject *)data)->len += sizeof(type); \
    } while (0)

static PyTypeObject BytesWriterType;

static bool
_grow_buffer(BytesWriterObject *data, Py_ssize_t n) {
    if (unlikely(n > PY_SSIZE_T_MAX - data->len)) {
        PyErr_NoMemory();
        return false;
    }
    Py_ssize_t target = data->len + n;
    Py_ssize_t size = data->capacity;
    do {
        if (unlikely(size > PY_SSIZE_T_MAX / 2)) {
            PyErr_NoMemory();
            return false;
        }
        size *= 2;
    } while (target >= size);
    char *new_buf;
    if (data->buf == data->data) {
        // Move from embedded buffer to heap-allocated buffer
        new_buf = PyMem_Malloc(size);
        if (new_buf != NULL) {
            memcpy(new_buf, data->data, data->len);
        }
    } else {
        new_buf = PyMem_Realloc(data->buf, size);
    }
    if (unlikely(new_buf == NULL)) {
        PyErr_NoMemory();
        return false;
    }
    data->buf = new_buf;
    data->capacity = size;
    return true;
}

static inline bool
ensure_bytes_writer_size(BytesWriterObject *data, Py_ssize_t n) {
    if (likely(data->capacity - data->len >= n)) {
        return true;
    } else {
        return _grow_buffer(data, n);
    }
}

static inline void
BytesWriter_init_internal(BytesWriterObject *self) {
    self->buf = self->data;
    self->len = 0;
    self->capacity = WRITER_EMBEDDED_BUF_LEN;
}

static PyObject*
BytesWriter_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    if (type != &BytesWriterType) {
        PyErr_SetString(PyExc_TypeError, "BytesWriter cannot be subclassed");
        return NULL;
    }

    BytesWriterObject *self = (BytesWriterObject *)type->tp_alloc(type, 0);
    if (self != NULL)
        BytesWriter_init_internal(self);
    return (PyObject *)self;
}

static PyObject *
BytesWriter_internal(void) {
    BytesWriterObject *self = (BytesWriterObject *)BytesWriterType.tp_alloc(&BytesWriterType, 0);
    if (self == NULL)
        return NULL;
    BytesWriter_init_internal(self);
    return (PyObject *)self;
}

static int
BytesWriter_init(BytesWriterObject *self, PyObject *args, PyObject *kwds)
{
    if (!PyArg_ParseTuple(args, "")) {
        return -1;
    }

    if (kwds != NULL && PyDict_Size(kwds) > 0) {
        PyErr_SetString(PyExc_TypeError,
                        "BytesWriter() takes no keyword arguments");
        return -1;
    }

    // Soft reset: free any heap buffer so re-initialization doesn't leak.
    if (self->buf != self->data && self->buf != NULL)
        PyMem_Free(self->buf);
    BytesWriter_init_internal(self);
    return 0;
}

static void
BytesWriter_dealloc(BytesWriterObject *self)
{
    if (self->buf != self->data) {
        PyMem_Free(self->buf);
        self->buf = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject*
BytesWriter_getvalue_internal(PyObject *self)
{
    BytesWriterObject *obj = (BytesWriterObject *)self;
    return PyBytes_FromStringAndSize(obj->buf, obj->len);
}

static PyObject*
BytesWriter_repr(BytesWriterObject *self)
{
    PyObject *value = BytesWriter_getvalue_internal((PyObject *)self);
    if (value == NULL) {
        return NULL;
    }
    PyObject *value_repr = PyObject_Repr(value);
    Py_DECREF(value);
    if (value_repr == NULL) {
        return NULL;
    }
    PyObject *result = PyUnicode_FromFormat("BytesWriter(%U)", value_repr);
    Py_DECREF(value_repr);
    return result;
}

static PyObject*
BytesWriter_getvalue(BytesWriterObject *self, PyObject *Py_UNUSED(ignored))
{
    return PyBytes_FromStringAndSize(self->buf, self->len);
}

static Py_ssize_t
BytesWriter_length(BytesWriterObject *self)
{
    return self->len;
}

static PyObject*
BytesWriter_item(BytesWriterObject *self, Py_ssize_t index)
{
    Py_ssize_t length = self->len;

    // Check bounds
    if (index < 0 || index >= length) {
        PyErr_SetString(PyExc_IndexError, "BytesWriter index out of range");
        return NULL;
    }

    // Return the byte at the given index as a Python int
    return PyLong_FromLong((unsigned char)self->buf[index]);
}

static int
BytesWriter_ass_item(BytesWriterObject *self, Py_ssize_t index, PyObject *value)
{
    Py_ssize_t length = self->len;

    // Check bounds
    if (index < 0 || index >= length) {
        PyErr_SetString(PyExc_IndexError, "BytesWriter index out of range");
        return -1;
    }

    // Check that value is not NULL (deletion not supported)
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "BytesWriter does not support item deletion");
        return -1;
    }

    // Convert value to uint8
    uint8_t byte_value = CPyLong_AsUInt8(value);
    if (unlikely(byte_value == CPY_LL_UINT_ERROR && PyErr_Occurred())) {
        CPy_TypeError("u8", value);
        return -1;
    }

    // Assign the byte
    self->buf[index] = (char)byte_value;
    return 0;
}

static PySequenceMethods BytesWriter_as_sequence = {
    .sq_length = (lenfunc)BytesWriter_length,
    .sq_item = (ssizeargfunc)BytesWriter_item,
    .sq_ass_item = (ssizeobjargproc)BytesWriter_ass_item,
};

static PyObject* BytesWriter_append(PyObject *self, PyObject *const *args, size_t nargs);
static PyObject* BytesWriter_write(PyObject *self, PyObject *const *args, size_t nargs);
static PyObject* BytesWriter_truncate(PyObject *self, PyObject *const *args, size_t nargs);

static PyMethodDef BytesWriter_methods[] = {
    {"append", (PyCFunction) BytesWriter_append, METH_FASTCALL,
     PyDoc_STR("Append a single byte to the buffer")
    },
    {"write", (PyCFunction) BytesWriter_write, METH_FASTCALL,
     PyDoc_STR("Append bytes to the buffer")
    },
    {"getvalue", (PyCFunction) BytesWriter_getvalue, METH_NOARGS,
     "Return the buffer content as bytes object"
    },
    {"truncate", (PyCFunction) BytesWriter_truncate, METH_FASTCALL,
     PyDoc_STR("Truncate the buffer to the specified size")
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject BytesWriterType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "BytesWriter",
    .tp_doc = PyDoc_STR("Memory buffer for building bytes objects from parts"),
    .tp_basicsize = sizeof(BytesWriterObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = BytesWriter_new,
    .tp_init = (initproc) BytesWriter_init,
    .tp_dealloc = (destructor) BytesWriter_dealloc,
    .tp_methods = BytesWriter_methods,
    .tp_as_sequence = &BytesWriter_as_sequence,
    .tp_repr = (reprfunc)BytesWriter_repr,
};

static inline bool
check_bytes_writer(PyObject *data) {
    if (unlikely(Py_TYPE(data) != &BytesWriterType)) {
        PyErr_Format(
            PyExc_TypeError, "data must be a BytesWriter object, got %s", Py_TYPE(data)->tp_name
        );
        return false;
    }
    return true;
}

static char
BytesWriter_write_internal(BytesWriterObject *self, PyObject *value) {
    const char *data;
    Py_ssize_t size;
    if (likely(PyBytes_Check(value))) {
        data = PyBytes_AS_STRING(value);
        size = PyBytes_GET_SIZE(value);
    } else {
        data = PyByteArray_AS_STRING(value);
        size = PyByteArray_GET_SIZE(value);
    }
    // Write bytes content.
    if (!ensure_bytes_writer_size(self, size))
        return CPY_NONE_ERROR;
    memcpy(self->buf + self->len, data, size);
    self->len += size;
    return CPY_NONE;
}

static PyObject*
BytesWriter_write(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 1)) {
        PyErr_Format(PyExc_TypeError,
                     "write() takes exactly 1 argument (%zu given)", nargs);
        return NULL;
    }
    if (!check_bytes_writer(self)) {
        return NULL;
    }
    PyObject *value = args[0];
    if (unlikely(!PyBytes_Check(value) && !PyByteArray_Check(value))) {
        PyErr_SetString(PyExc_TypeError, "value must be a bytes or bytearray object");
        return NULL;
    }
    if (unlikely(BytesWriter_write_internal((BytesWriterObject *)self, value) == CPY_NONE_ERROR)) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static inline char
BytesWriter_append_internal(BytesWriterObject *self, uint8_t value) {
    if (!ensure_bytes_writer_size(self, 1))
        return CPY_NONE_ERROR;
    _WRITE_BYTES(self, uint8_t, value);
    return CPY_NONE;
}

static PyObject*
BytesWriter_append(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 1)) {
        PyErr_Format(PyExc_TypeError,
                     "append() takes exactly 1 argument (%zu given)", nargs);
        return NULL;
    }
    if (!check_bytes_writer(self)) {
        return NULL;
    }
    PyObject *value = args[0];
    uint8_t unboxed = CPyLong_AsUInt8(value);
    if (unlikely(unboxed == CPY_LL_UINT_ERROR && PyErr_Occurred())) {
        CPy_TypeError("u8", value);
        return NULL;
    }
    if (unlikely(BytesWriter_append_internal((BytesWriterObject *)self, unboxed) == CPY_NONE_ERROR)) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static char
BytesWriter_truncate_internal(PyObject *self, int64_t size) {
    BytesWriterObject *writer = (BytesWriterObject *)self;
    Py_ssize_t current_size = writer->len;

    // Validate size is non-negative
    if (size < 0) {
        PyErr_SetString(PyExc_ValueError, "size must be non-negative");
        return CPY_NONE_ERROR;
    }

    // Validate size doesn't exceed current size
    if (size > current_size) {
        PyErr_SetString(PyExc_ValueError, "size cannot be larger than current buffer size");
        return CPY_NONE_ERROR;
    }

    writer->len = size;
    return CPY_NONE;
}

static PyObject*
BytesWriter_truncate(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 1)) {
        PyErr_Format(PyExc_TypeError,
                     "truncate() takes exactly 1 argument (%zu given)", nargs);
        return NULL;
    }
    if (!check_bytes_writer(self)) {
        return NULL;
    }

    PyObject *size_obj = args[0];
    int overflow;
    long long size = PyLong_AsLongLongAndOverflow(size_obj, &overflow);

    if (size == -1 && PyErr_Occurred()) {
        return NULL;
    }
    if (overflow != 0) {
        PyErr_SetString(PyExc_ValueError, "integer out of range");
        return NULL;
    }

    if (unlikely(BytesWriter_truncate_internal(self, size) == CPY_NONE_ERROR)) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyTypeObject *
BytesWriter_type_internal(void) {
    return &BytesWriterType;  // Return borrowed reference
};

static CPyTagged
BytesWriter_len_internal(PyObject *self) {
    BytesWriterObject *writer = (BytesWriterObject *)self;
    return writer->len << 1;
}

//
// StringWriter
//

static PyTypeObject StringWriterType;

static void convert_string_data_in_place(char *buf, Py_ssize_t len,
                                         char old_kind, char new_kind);

// Helper to grow string buffer and optionally convert to new kind
// Returns true on success, false on failure (with PyErr set)
// Updates self->buf, self->capacity, and self->kind
static bool
grow_string_buffer_helper(StringWriterObject *self, Py_ssize_t target_capacity, char new_kind) {
    char old_kind = self->kind;
    Py_ssize_t new_capacity = self->capacity;

    // Limit so that (new_capacity * 2) * new_kind stays within Py_ssize_t.
    // new_kind is always 1, 2, or 4, so use a shift instead of division.
    int shift = (new_kind == 4) ? 3 : (new_kind == 2 ? 2 : 1);
    Py_ssize_t cap_limit = PY_SSIZE_T_MAX >> shift;
    while (target_capacity >= new_capacity) {
        if (unlikely(new_capacity > cap_limit)) {
            PyErr_NoMemory();
            return false;
        }
        new_capacity *= 2;
    }

    if (unlikely(new_capacity > cap_limit * 2)) {
        PyErr_NoMemory();
        return false;
    }

    Py_ssize_t size_bytes = new_capacity * new_kind;
    char *new_buf;
    bool from_embedded = (self->buf == self->data);

    if (from_embedded) {
        // Move from embedded buffer to heap-allocated buffer
        new_buf = PyMem_Malloc(size_bytes);
        if (new_buf != NULL) {
            // Copy existing data from embedded buffer
            memcpy(new_buf, self->data, self->len * old_kind);
        }
    } else {
        // Realloc existing heap buffer
        new_buf = PyMem_Realloc(self->buf, size_bytes);
    }

    if (unlikely(new_buf == NULL)) {
        PyErr_NoMemory();
        return false;
    }

    // Convert data if kind changed
    if (old_kind != new_kind) {
        convert_string_data_in_place(new_buf, self->len, old_kind, new_kind);
    }

    self->buf = new_buf;
    self->capacity = new_capacity;
    self->kind = new_kind;
    return true;
}

static bool grow_string_buffer(StringWriterObject *data, Py_ssize_t n) {
    if (unlikely(n > PY_SSIZE_T_MAX - data->len)) {
        PyErr_NoMemory();
        return false;
    }
    return grow_string_buffer_helper(data, data->len + n, data->kind);
}

static inline bool
ensure_string_writer_size(StringWriterObject *data, Py_ssize_t n) {
    if (likely(data->capacity - data->len >= n)) {
        return true;
    } else {
        // Don't inline the grow function since this is slow path and we
        // want to keep this as short as possible for better inlining
        return grow_string_buffer(data, n);
    }
}

static inline void
StringWriter_init_internal(StringWriterObject *self) {
    self->buf = self->data;
    self->kind = 1;
    self->len = 0;
    self->capacity = WRITER_EMBEDDED_BUF_LEN;
}

static PyObject*
StringWriter_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    if (type != &StringWriterType) {
        PyErr_SetString(PyExc_TypeError, "StringWriter cannot be subclassed");
        return NULL;
    }

    StringWriterObject *self = (StringWriterObject *)type->tp_alloc(type, 0);
    if (self != NULL)
        StringWriter_init_internal(self);
    return (PyObject *)self;
}

static PyObject *
StringWriter_internal(void) {
    StringWriterObject *self = (StringWriterObject *)StringWriterType.tp_alloc(&StringWriterType, 0);
    if (self == NULL)
        return NULL;
    StringWriter_init_internal(self);
    return (PyObject *)self;
}

static int
StringWriter_init(StringWriterObject *self, PyObject *args, PyObject *kwds)
{
    if (!PyArg_ParseTuple(args, "")) {
        return -1;
    }

    if (kwds != NULL && PyDict_Size(kwds) > 0) {
        PyErr_SetString(PyExc_TypeError,
                        "StringWriter() takes no keyword arguments");
        return -1;
    }

    // Soft reset: free any heap buffer so re-initialization doesn't leak.
    if (self->buf != self->data && self->buf != NULL)
        PyMem_Free(self->buf);
    StringWriter_init_internal(self);
    return 0;
}

static void
StringWriter_dealloc(StringWriterObject *self)
{
    if (self->buf != self->data) {
        PyMem_Free(self->buf);
        self->buf = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject*
StringWriter_getvalue_internal(PyObject *self)
{
    StringWriterObject *obj = (StringWriterObject *)self;
    return PyUnicode_FromKindAndData(obj->kind, obj->buf, obj->len);
}

static PyObject*
StringWriter_repr(StringWriterObject *self)
{
    PyObject *value = StringWriter_getvalue_internal((PyObject *)self);
    if (value == NULL) {
        return NULL;
    }
    PyObject *value_repr = PyObject_Repr(value);
    Py_DECREF(value);
    if (value_repr == NULL) {
        return NULL;
    }
    PyObject *result = PyUnicode_FromFormat("StringWriter(%U)", value_repr);
    Py_DECREF(value_repr);
    return result;
}

static PyObject*
StringWriter_getvalue(StringWriterObject *self, PyObject *Py_UNUSED(ignored))
{
    return PyUnicode_FromKindAndData(self->kind, self->buf, self->len);
}

static Py_ssize_t
StringWriter_length(StringWriterObject *self)
{
    return self->len;
}

static PyObject*
StringWriter_item(StringWriterObject *self, Py_ssize_t index)
{
    Py_ssize_t length = self->len;

    // Check bounds
    if (index < 0 || index >= length) {
        PyErr_SetString(PyExc_IndexError, "StringWriter index out of range");
        return NULL;
    }

    // Read the character at the given index based on kind using memcpy
    uint32_t value;
    if (self->kind == 1) {
        uint8_t val;
        memcpy(&val, self->buf + index, 1);
        value = val;
    } else if (self->kind == 2) {
        uint16_t val;
        memcpy(&val, self->buf + index * 2, 2);
        value = val;
    } else {
        memcpy(&value, self->buf + index * 4, 4);
    }
    return PyLong_FromLong(value);
}

static PySequenceMethods StringWriter_as_sequence = {
    .sq_length = (lenfunc)StringWriter_length,
    .sq_item = (ssizeargfunc)StringWriter_item,
};

static PyObject* StringWriter_append(PyObject *self, PyObject *const *args, size_t nargs);
static PyObject* StringWriter_write(PyObject *self, PyObject *const *args, size_t nargs);

static PyMethodDef StringWriter_methods[] = {
    {"append", (PyCFunction) StringWriter_append, METH_FASTCALL,
     PyDoc_STR("Append a single character (as int codepoint) to the buffer")
    },
    {"write", (PyCFunction) StringWriter_write, METH_FASTCALL,
     PyDoc_STR("Append a string to the buffer")
    },
    {"getvalue", (PyCFunction) StringWriter_getvalue, METH_NOARGS,
     "Return the buffer content as str object"
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject StringWriterType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "StringWriter",
    .tp_doc = PyDoc_STR("Memory buffer for building string objects from parts"),
    .tp_basicsize = sizeof(StringWriterObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = StringWriter_new,
    .tp_init = (initproc) StringWriter_init,
    .tp_dealloc = (destructor) StringWriter_dealloc,
    .tp_methods = StringWriter_methods,
    .tp_as_sequence = &StringWriter_as_sequence,
    .tp_repr = (reprfunc)StringWriter_repr,
};

static inline bool
check_string_writer(PyObject *data) {
    if (unlikely(Py_TYPE(data) != &StringWriterType)) {
        PyErr_Format(
            PyExc_TypeError, "data must be a StringWriter object, got %s", Py_TYPE(data)->tp_name
        );
        return false;
    }
    return true;
}

static char string_writer_switch_kind(StringWriterObject *self, int32_t value);

static char
StringWriter_write_internal(PyObject *obj, PyObject *value) {
    StringWriterObject *self = (StringWriterObject *)obj;
    Py_ssize_t str_len = PyUnicode_GET_LENGTH(value);
    if (str_len == 0) {
        return CPY_NONE;
    }

    int src_kind = PyUnicode_KIND(value);
    void *src_data = PyUnicode_DATA(value);
    int self_kind = self->kind;

    // Switch kind if source requires wider characters
    if (src_kind > self_kind) {
        // Use value in the source kind range to trigger proper kind switch
        int32_t codepoint = (src_kind == 2) ? 0x100 : 0x10000;
        if (string_writer_switch_kind(self, codepoint) == CPY_NONE_ERROR) {
            return CPY_NONE_ERROR;
        }
        self_kind = self->kind;
    }

    // Ensure we have enough space
    if (!ensure_string_writer_size(self, str_len)) {
        return CPY_NONE_ERROR;
    }

    // Copy data - ASCII/Latin1 (kind 1) are handled uniformly
    if (self_kind == src_kind) {
        // Same kind, direct copy
        memcpy(self->buf + self->len * self_kind, src_data, str_len * src_kind);
    } else {
        // Different kinds, convert character by character
        for (Py_ssize_t i = 0; i < str_len; i++) {
            Py_UCS4 ch = PyUnicode_READ(src_kind, src_data, i);
            PyUnicode_WRITE(self_kind, self->buf, self->len + i, ch);
        }
    }

    self->len += str_len;
    return CPY_NONE;
}

static PyObject*
StringWriter_write(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 1)) {
        PyErr_Format(PyExc_TypeError,
                     "write() takes exactly 1 argument (%zu given)", nargs);
        return NULL;
    }
    if (!check_string_writer(self)) {
        return NULL;
    }
    PyObject *value = args[0];
    if (unlikely(!PyUnicode_Check(value))) {
        PyErr_SetString(PyExc_TypeError, "value must be a str object");
        return NULL;
    }
    if (unlikely(StringWriter_write_internal(self, value) == CPY_NONE_ERROR)) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// Convert string data to next larger kind (1->2 or 2->4)
static void convert_string_data_in_place(char *buf, Py_ssize_t len,
                                         char old_kind, char new_kind) {
    if (old_kind == 1 && new_kind == 2) {
        // Convert backwards to avoid overwriting
        for (Py_ssize_t i = len - 1; i >= 0; i--) {
            uint8_t val = (uint8_t)buf[i];
            uint16_t expanded = val;
            memcpy(buf + i * 2, &expanded, 2);
        }
    } else if (old_kind == 2 && new_kind == 4) {
        // Convert backwards to avoid overwriting
        for (Py_ssize_t i = len - 1; i >= 0; i--) {
            uint16_t val;
            memcpy(&val, buf + i * 2, 2);
            uint32_t expanded = val;
            memcpy(buf + i * 4, &expanded, 4);
        }
    } else {
        assert(false);
    }
}

static char convert_string_buffer_kind(StringWriterObject *self, char old_kind, char new_kind) {
    // new_kind is always 2 or 4, so the max representable len is PY_SSIZE_T_MAX >> {1,2}.
    // Use a shift instead of a division for cheap overflow check.
    Py_ssize_t max_len = PY_SSIZE_T_MAX >> (new_kind == 4 ? 2 : 1);
    if (unlikely(self->len > max_len)) {
        PyErr_NoMemory();
        return CPY_NONE_ERROR;
    }
    // Current buffer size in bytes
    Py_ssize_t current_buf_size = (self->buf == self->data) ? WRITER_EMBEDDED_BUF_LEN : (self->capacity * old_kind);
    // Needed buffer size in bytes for new kind
    Py_ssize_t needed_size = self->len * new_kind;

    if (current_buf_size >= needed_size) {
        // Convert in place
        convert_string_data_in_place(self->buf, self->len, old_kind, new_kind);
        self->kind = new_kind;
        self->capacity = current_buf_size / new_kind;
    } else {
        // Need to allocate new buffer
        if (!grow_string_buffer_helper(self, self->len, new_kind)) {
            return CPY_NONE_ERROR;
        }
    }
    return CPY_NONE;
}

static char string_writer_switch_kind(StringWriterObject *self, int32_t value) {
    if (self->kind == 1) {
        // Either kind 1 -> 2 or 1 -> 4. First switch to kind 2.
        if (convert_string_buffer_kind(self, 1, 2) == CPY_NONE_ERROR)
            return CPY_NONE_ERROR;
        if ((uint32_t)value > 0xffff) {
            // Call recursively to switch from kind 2 to 4
            return string_writer_switch_kind(self, value);
        }
        return CPY_NONE;
    } else {
        // Must be kind 2 -> 4
        assert(self->kind == 2);
        assert((uint32_t)value > 0xffff);
        return convert_string_buffer_kind(self, 2, 4);
    }
}

// Handle all append cases except for append that stays within kind 1
static char string_append_slow_path(StringWriterObject *self, int32_t value) {
    if (self->kind == 2) {
        if ((uint32_t)value <= 0xffff) {
            // Fast path - kind 2 stays the same
            if (!ensure_string_writer_size(self, 1))
                return CPY_NONE_ERROR;
            // Copy 2-byte character to buffer
            uint16_t val16 = (uint16_t)value;
            memcpy(self->buf + self->len * 2, &val16, 2);
            self->len++;
            return CPY_NONE;
        }
        // Always validate code point range before promotion (but after fast path).
        if (unlikely((uint32_t)value > 0x10FFFF))
            goto fail_range;
        if (string_writer_switch_kind(self, value) == CPY_NONE_ERROR)
            return CPY_NONE_ERROR;
        return string_append_slow_path(self, value);
    } else if (self->kind == 1) {
        // Check precondition -- this must only be used on slow path
        assert((uint32_t)value > 0xff);
        if (unlikely((uint32_t)value > 0x10FFFF))
            goto fail_range;
        if (string_writer_switch_kind(self, value) == CPY_NONE_ERROR)
            return CPY_NONE_ERROR;
        return string_append_slow_path(self, value);
    }
    assert(self->kind == 4);
    if (unlikely((uint32_t)value > 0x10FFFF))
        goto fail_range;
    if (!ensure_string_writer_size(self, 1))
        return CPY_NONE_ERROR;
    // Copy 4-byte character to buffer
    uint32_t val32 = (uint32_t)value;
    memcpy(self->buf + self->len * 4, &val32, 4);
    self->len++;
    return CPY_NONE;

fail_range:
    PyErr_Format(PyExc_ValueError,
                "code point %d is outside valid Unicode range (0-1114111)", value);
    return CPY_NONE_ERROR;
}

static inline char
StringWriter_append_internal(StringWriterObject *self, int32_t value) {
    char kind = self->kind;
    if (kind == 1 && (uint32_t)value < 256) {
        if (!ensure_string_writer_size(self, 1))
            return CPY_NONE_ERROR;
        self->buf[self->len++] = value;
        self->kind = kind;
        return CPY_NONE;
    }
    return string_append_slow_path(self, value);
}

static PyObject*
StringWriter_append(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 1)) {
        PyErr_Format(PyExc_TypeError,
                     "append() takes exactly 1 argument (%zu given)", nargs);
        return NULL;
    }
    if (!check_string_writer(self)) {
        return NULL;
    }
    PyObject *value = args[0];
    int32_t unboxed = CPyLong_AsInt32(value);
    if (unlikely(unboxed == CPY_LL_INT_ERROR && PyErr_Occurred())) {
        CPy_TypeError("i32", value);
        return NULL;
    }
    if (unlikely(StringWriter_append_internal((StringWriterObject *)self, unboxed) == CPY_NONE_ERROR)) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyTypeObject *
StringWriter_type_internal(void) {
    return &StringWriterType;  // Return borrowed reference
};

static CPyTagged
StringWriter_len_internal(PyObject *self) {
    StringWriterObject *writer = (StringWriterObject *)self;
    return writer->len << 1;
}

// End of StringWriter

// Helper for write_*_le/be functions - validates args and returns BytesWriter
static inline BytesWriterObject *
parse_write_args(PyObject *const *args, size_t nargs, const char *func_name) {
    if (unlikely(nargs != 2)) {
        PyErr_Format(PyExc_TypeError,
                     "%s() takes exactly 2 arguments (%zu given)", func_name, nargs);
        return NULL;
    }
    PyObject *writer = args[0];
    if (!check_bytes_writer(writer)) {
        return NULL;
    }
    return (BytesWriterObject *)writer;
}

// Helper for read_*_le/be functions - validates args and returns data pointer
// Returns NULL on error, sets *out_index to the validated index on success
static inline const unsigned char *
parse_read_args(PyObject *const *args, size_t nargs, const char *func_name,
                    Py_ssize_t num_bytes, int64_t *out_index) {
    if (unlikely(nargs != 2)) {
        PyErr_Format(PyExc_TypeError,
                     "%s() takes exactly 2 arguments (%zu given)", func_name, nargs);
        return NULL;
    }
    PyObject *bytes_obj = args[0];
    if (unlikely(!PyBytes_Check(bytes_obj))) {
        PyErr_Format(PyExc_TypeError, "%s() argument 1 must be bytes", func_name);
        return NULL;
    }
    int64_t index = CPyLong_AsInt64(args[1]);
    if (unlikely(index == CPY_LL_INT_ERROR && PyErr_Occurred())) {
        return NULL;
    }
    if (unlikely(index < 0)) {
        PyErr_SetString(PyExc_ValueError, "index must be non-negative");
        return NULL;
    }
    Py_ssize_t size = PyBytes_GET_SIZE(bytes_obj);
    if (unlikely(index > size - num_bytes)) {
        PyErr_Format(PyExc_IndexError,
                     "index %lld out of range for bytes of length %zd",
                     (long long)index, size);
        return NULL;
    }
    *out_index = index;
    return (const unsigned char *)PyBytes_AS_STRING(bytes_obj);
}

static PyObject*
write_i16_le(PyObject *module, PyObject *const *args, size_t nargs) {
    BytesWriterObject *bw = parse_write_args(args, nargs, "write_i16_le");
    if (bw == NULL)
        return NULL;
    int16_t unboxed = CPyLong_AsInt16(args[1]);
    if (unlikely(unboxed == CPY_LL_INT_ERROR && PyErr_Occurred()))
        return NULL;
    if (unlikely(!ensure_bytes_writer_size(bw, 2)))
        return NULL;
    BytesWriter_WriteI16LEUnsafe(bw, unboxed);
    Py_RETURN_NONE;
}

static PyObject*
write_i16_be(PyObject *module, PyObject *const *args, size_t nargs) {
    BytesWriterObject *bw = parse_write_args(args, nargs, "write_i16_be");
    if (bw == NULL)
        return NULL;
    int16_t unboxed = CPyLong_AsInt16(args[1]);
    if (unlikely(unboxed == CPY_LL_INT_ERROR && PyErr_Occurred()))
        return NULL;
    if (unlikely(!ensure_bytes_writer_size(bw, 2)))
        return NULL;
    BytesWriter_WriteI16BEUnsafe(bw, unboxed);
    Py_RETURN_NONE;
}

static PyObject*
read_i16_le(PyObject *module, PyObject *const *args, size_t nargs) {
    int64_t index;
    const unsigned char *data = parse_read_args(args, nargs, "read_i16_le", 2, &index);
    if (data == NULL)
        return NULL;
    return PyLong_FromLong(CPyBytes_ReadI16LEUnsafe(data + index));
}

static PyObject*
read_i16_be(PyObject *module, PyObject *const *args, size_t nargs) {
    int64_t index;
    const unsigned char *data = parse_read_args(args, nargs, "read_i16_be", 2, &index);
    if (data == NULL)
        return NULL;
    return PyLong_FromLong(CPyBytes_ReadI16BEUnsafe(data + index));
}

static PyObject*
write_i32_le(PyObject *module, PyObject *const *args, size_t nargs) {
    BytesWriterObject *bw = parse_write_args(args, nargs, "write_i32_le");
    if (bw == NULL)
        return NULL;
    int32_t unboxed = CPyLong_AsInt32(args[1]);
    if (unlikely(unboxed == CPY_LL_INT_ERROR && PyErr_Occurred()))
        return NULL;
    if (unlikely(!ensure_bytes_writer_size(bw, 4)))
        return NULL;
    BytesWriter_WriteI32LEUnsafe(bw, unboxed);
    Py_RETURN_NONE;
}

static PyObject*
write_i32_be(PyObject *module, PyObject *const *args, size_t nargs) {
    BytesWriterObject *bw = parse_write_args(args, nargs, "write_i32_be");
    if (bw == NULL)
        return NULL;
    int32_t unboxed = CPyLong_AsInt32(args[1]);
    if (unlikely(unboxed == CPY_LL_INT_ERROR && PyErr_Occurred()))
        return NULL;
    if (unlikely(!ensure_bytes_writer_size(bw, 4)))
        return NULL;
    BytesWriter_WriteI32BEUnsafe(bw, unboxed);
    Py_RETURN_NONE;
}

static PyObject*
read_i32_le(PyObject *module, PyObject *const *args, size_t nargs) {
    int64_t index;
    const unsigned char *data = parse_read_args(args, nargs, "read_i32_le", 4, &index);
    if (data == NULL)
        return NULL;
    return PyLong_FromLong(CPyBytes_ReadI32LEUnsafe(data + index));
}

static PyObject*
read_i32_be(PyObject *module, PyObject *const *args, size_t nargs) {
    int64_t index;
    const unsigned char *data = parse_read_args(args, nargs, "read_i32_be", 4, &index);
    if (data == NULL)
        return NULL;
    return PyLong_FromLong(CPyBytes_ReadI32BEUnsafe(data + index));
}

static PyObject*
write_i64_le(PyObject *module, PyObject *const *args, size_t nargs) {
    BytesWriterObject *bw = parse_write_args(args, nargs, "write_i64_le");
    if (bw == NULL)
        return NULL;
    int64_t unboxed = CPyLong_AsInt64(args[1]);
    if (unlikely(unboxed == CPY_LL_INT_ERROR && PyErr_Occurred()))
        return NULL;
    if (unlikely(!ensure_bytes_writer_size(bw, 8)))
        return NULL;
    BytesWriter_WriteI64LEUnsafe(bw, unboxed);
    Py_RETURN_NONE;
}

static PyObject*
write_i64_be(PyObject *module, PyObject *const *args, size_t nargs) {
    BytesWriterObject *bw = parse_write_args(args, nargs, "write_i64_be");
    if (bw == NULL)
        return NULL;
    int64_t unboxed = CPyLong_AsInt64(args[1]);
    if (unlikely(unboxed == CPY_LL_INT_ERROR && PyErr_Occurred()))
        return NULL;
    if (unlikely(!ensure_bytes_writer_size(bw, 8)))
        return NULL;
    BytesWriter_WriteI64BEUnsafe(bw, unboxed);
    Py_RETURN_NONE;
}

static PyObject*
read_i64_le(PyObject *module, PyObject *const *args, size_t nargs) {
    int64_t index;
    const unsigned char *data = parse_read_args(args, nargs, "read_i64_le", 8, &index);
    if (data == NULL)
        return NULL;
    return PyLong_FromLongLong(CPyBytes_ReadI64LEUnsafe(data + index));
}

static PyObject*
read_i64_be(PyObject *module, PyObject *const *args, size_t nargs) {
    int64_t index;
    const unsigned char *data = parse_read_args(args, nargs, "read_i64_be", 8, &index);
    if (data == NULL)
        return NULL;
    return PyLong_FromLongLong(CPyBytes_ReadI64BEUnsafe(data + index));
}

static PyObject*
write_f32_le(PyObject *module, PyObject *const *args, size_t nargs) {
    BytesWriterObject *bw = parse_write_args(args, nargs, "write_f32_le");
    if (bw == NULL)
        return NULL;
    double unboxed = PyFloat_AsDouble(args[1]);
    if (unlikely(unboxed == -1.0 && PyErr_Occurred()))
        return NULL;
    if (unlikely(!ensure_bytes_writer_size(bw, 4)))
        return NULL;
    BytesWriter_WriteF32LEUnsafe(bw, (float)unboxed);
    Py_RETURN_NONE;
}

static PyObject*
write_f32_be(PyObject *module, PyObject *const *args, size_t nargs) {
    BytesWriterObject *bw = parse_write_args(args, nargs, "write_f32_be");
    if (bw == NULL)
        return NULL;
    double unboxed = PyFloat_AsDouble(args[1]);
    if (unlikely(unboxed == -1.0 && PyErr_Occurred()))
        return NULL;
    if (unlikely(!ensure_bytes_writer_size(bw, 4)))
        return NULL;
    BytesWriter_WriteF32BEUnsafe(bw, (float)unboxed);
    Py_RETURN_NONE;
}

static PyObject*
read_f32_le(PyObject *module, PyObject *const *args, size_t nargs) {
    int64_t index;
    const unsigned char *data = parse_read_args(args, nargs, "read_f32_le", 4, &index);
    if (data == NULL)
        return NULL;
    return PyFloat_FromDouble((double)CPyBytes_ReadF32LEUnsafe(data + index));
}

static PyObject*
read_f32_be(PyObject *module, PyObject *const *args, size_t nargs) {
    int64_t index;
    const unsigned char *data = parse_read_args(args, nargs, "read_f32_be", 4, &index);
    if (data == NULL)
        return NULL;
    return PyFloat_FromDouble((double)CPyBytes_ReadF32BEUnsafe(data + index));
}

static PyObject*
write_f64_le(PyObject *module, PyObject *const *args, size_t nargs) {
    BytesWriterObject *bw = parse_write_args(args, nargs, "write_f64_le");
    if (bw == NULL)
        return NULL;
    double unboxed = PyFloat_AsDouble(args[1]);
    if (unlikely(unboxed == -1.0 && PyErr_Occurred()))
        return NULL;
    if (unlikely(!ensure_bytes_writer_size(bw, 8)))
        return NULL;
    BytesWriter_WriteF64LEUnsafe(bw, unboxed);
    Py_RETURN_NONE;
}

static PyObject*
write_f64_be(PyObject *module, PyObject *const *args, size_t nargs) {
    BytesWriterObject *bw = parse_write_args(args, nargs, "write_f64_be");
    if (bw == NULL)
        return NULL;
    double unboxed = PyFloat_AsDouble(args[1]);
    if (unlikely(unboxed == -1.0 && PyErr_Occurred()))
        return NULL;
    if (unlikely(!ensure_bytes_writer_size(bw, 8)))
        return NULL;
    BytesWriter_WriteF64BEUnsafe(bw, unboxed);
    Py_RETURN_NONE;
}

static PyObject*
read_f64_le(PyObject *module, PyObject *const *args, size_t nargs) {
    int64_t index;
    const unsigned char *data = parse_read_args(args, nargs, "read_f64_le", 8, &index);
    if (data == NULL)
        return NULL;
    return PyFloat_FromDouble(CPyBytes_ReadF64LEUnsafe(data + index));
}

static PyObject*
read_f64_be(PyObject *module, PyObject *const *args, size_t nargs) {
    int64_t index;
    const unsigned char *data = parse_read_args(args, nargs, "read_f64_be", 8, &index);
    if (data == NULL)
        return NULL;
    return PyFloat_FromDouble(CPyBytes_ReadF64BEUnsafe(data + index));
}

static PyMethodDef librt_strings_module_methods[] = {
    {"write_i16_le", (PyCFunction) write_i16_le, METH_FASTCALL,
     PyDoc_STR("Write a 16-bit signed integer to BytesWriter in little-endian format")
    },
    {"write_i16_be", (PyCFunction) write_i16_be, METH_FASTCALL,
     PyDoc_STR("Write a 16-bit signed integer to BytesWriter in big-endian format")
    },
    {"read_i16_le", (PyCFunction) read_i16_le, METH_FASTCALL,
     PyDoc_STR("Read a 16-bit signed integer from bytes in little-endian format")
    },
    {"read_i16_be", (PyCFunction) read_i16_be, METH_FASTCALL,
     PyDoc_STR("Read a 16-bit signed integer from bytes in big-endian format")
    },
    {"write_i32_le", (PyCFunction) write_i32_le, METH_FASTCALL,
     PyDoc_STR("Write a 32-bit signed integer to BytesWriter in little-endian format")
    },
    {"write_i32_be", (PyCFunction) write_i32_be, METH_FASTCALL,
     PyDoc_STR("Write a 32-bit signed integer to BytesWriter in big-endian format")
    },
    {"read_i32_le", (PyCFunction) read_i32_le, METH_FASTCALL,
     PyDoc_STR("Read a 32-bit signed integer from bytes in little-endian format")
    },
    {"read_i32_be", (PyCFunction) read_i32_be, METH_FASTCALL,
     PyDoc_STR("Read a 32-bit signed integer from bytes in big-endian format")
    },
    {"write_i64_le", (PyCFunction) write_i64_le, METH_FASTCALL,
     PyDoc_STR("Write a 64-bit signed integer to BytesWriter in little-endian format")
    },
    {"write_i64_be", (PyCFunction) write_i64_be, METH_FASTCALL,
     PyDoc_STR("Write a 64-bit signed integer to BytesWriter in big-endian format")
    },
    {"read_i64_le", (PyCFunction) read_i64_le, METH_FASTCALL,
     PyDoc_STR("Read a 64-bit signed integer from bytes in little-endian format")
    },
    {"read_i64_be", (PyCFunction) read_i64_be, METH_FASTCALL,
     PyDoc_STR("Read a 64-bit signed integer from bytes in big-endian format")
    },
    {"write_f32_le", (PyCFunction) write_f32_le, METH_FASTCALL,
     PyDoc_STR("Write a 32-bit float to BytesWriter in little-endian format")
    },
    {"write_f32_be", (PyCFunction) write_f32_be, METH_FASTCALL,
     PyDoc_STR("Write a 32-bit float to BytesWriter in big-endian format")
    },
    {"read_f32_le", (PyCFunction) read_f32_le, METH_FASTCALL,
     PyDoc_STR("Read a 32-bit float from bytes in little-endian format")
    },
    {"read_f32_be", (PyCFunction) read_f32_be, METH_FASTCALL,
     PyDoc_STR("Read a 32-bit float from bytes in big-endian format")
    },
    {"write_f64_le", (PyCFunction) write_f64_le, METH_FASTCALL,
     PyDoc_STR("Write a 64-bit float to BytesWriter in little-endian format")
    },
    {"write_f64_be", (PyCFunction) write_f64_be, METH_FASTCALL,
     PyDoc_STR("Write a 64-bit float to BytesWriter in big-endian format")
    },
    {"read_f64_le", (PyCFunction) read_f64_le, METH_FASTCALL,
     PyDoc_STR("Read a 64-bit float from bytes in little-endian format")
    },
    {"read_f64_be", (PyCFunction) read_f64_be, METH_FASTCALL,
     PyDoc_STR("Read a 64-bit float from bytes in big-endian format")
    },
    {NULL, NULL, 0, NULL}
};

static int
strings_abi_version(void) {
    return LIBRT_STRINGS_ABI_VERSION;
}

static int
strings_api_version(void) {
    return LIBRT_STRINGS_API_VERSION;
}

static int
librt_strings_module_exec(PyObject *m)
{
    if (PyType_Ready(&BytesWriterType) < 0) {
        return -1;
    }
    if (PyType_Ready(&StringWriterType) < 0) {
        return -1;
    }
    if (PyModule_AddObjectRef(m, "BytesWriter", (PyObject *) &BytesWriterType) < 0) {
        return -1;
    }
    if (PyModule_AddObjectRef(m, "StringWriter", (PyObject *) &StringWriterType) < 0) {
        return -1;
    }

    // Export mypy internal C API, be careful with the order!
    static void *librt_strings_api[LIBRT_STRINGS_API_LEN] = {
        (void *)strings_abi_version,
        (void *)strings_api_version,
        (void *)BytesWriter_internal,
        (void *)BytesWriter_getvalue_internal,
        (void *)BytesWriter_append_internal,
        (void *)_grow_buffer,
        (void *)BytesWriter_type_internal,
        (void *)BytesWriter_truncate_internal,
        (void *)StringWriter_internal,
        (void *)StringWriter_getvalue_internal,
        (void *)string_append_slow_path,
        (void *)StringWriter_type_internal,
        (void *)StringWriter_write_internal,
        (void *)grow_string_buffer,
    };
    PyObject *c_api_object = PyCapsule_New((void *)librt_strings_api, "librt.strings._C_API", NULL);
    if (PyModule_Add(m, "_C_API", c_api_object) < 0) {
        return -1;
    }
    return 0;
}

static PyModuleDef_Slot librt_strings_module_slots[] = {
    {Py_mod_exec, librt_strings_module_exec},
#ifdef Py_MOD_GIL_NOT_USED
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};

static PyModuleDef librt_strings_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "strings",
    .m_doc = "Utilities for working with str and bytes objects",
    .m_size = 0,
    .m_methods = librt_strings_module_methods,
    .m_slots = librt_strings_module_slots,
};

PyMODINIT_FUNC
PyInit_strings(void)
{
    intern_strings();
    return PyModuleDef_Init(&librt_strings_module);
}
