#include "pythoncapi_compat.h"

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdint.h>
#include <string.h>
#include "CPy.h"
#define LIBRT_INTERNAL_MODULE
#include "librt_internal.h"

#define START_SIZE 512

// See comment in read_int_internal() on motivation for these values.
#define MIN_ONE_BYTE_INT -10
#define MAX_ONE_BYTE_INT 117  // 2 ** 7 - 1 - 10
#define MIN_TWO_BYTES_INT -100
#define MAX_TWO_BYTES_INT 16283  // 2 ** (8 + 6) - 1 - 100
#define MIN_FOUR_BYTES_INT -10000
#define MAX_FOUR_BYTES_INT 536860911  // 2 ** (3 * 8 + 5) - 1 - 10000

#define TWO_BYTES_INT_BIT 1
#define FOUR_BYTES_INT_BIT 2
#define LONG_INT_BIT 4

#define FOUR_BYTES_INT_TRAILER 3
// We add one reserved bit here so that we can potentially support
// 8 bytes format in the future.
#define LONG_INT_TRAILER 15

#define CPY_BOOL_ERROR 2

#define _CHECK_READ_BUFFER(data, err)  if (unlikely(_check_read_buffer(data) == CPY_NONE_ERROR)) \
                                           return err;
#define _CHECK_WRITE_BUFFER(data, err) if (unlikely(_check_write_buffer(data) == CPY_NONE_ERROR)) \
                                           return err;
#define _CHECK_WRITE(data, need)        if (unlikely(_check_size((WriteBufferObject *)data, need) == CPY_NONE_ERROR)) \
                                           return CPY_NONE_ERROR;
#define _CHECK_READ(data, size, err)   if (unlikely(_check_read((ReadBufferObject *)data, size) == CPY_NONE_ERROR)) \
                                           return err;

#define _READ(result, data, type) \
    do { \
        memcpy((void *) result, ((ReadBufferObject *)data)->ptr, sizeof(type)); \
        ((ReadBufferObject *)data)->ptr += sizeof(type); \
    } while (0)

#define _WRITE(data, type, v) \
    do { \
       type temp = v; \
       memcpy(((WriteBufferObject *)data)->ptr, (const void *) &temp, sizeof(type)); \
       ((WriteBufferObject *)data)->ptr += sizeof(type); \
    } while (0)

//
// ReadBuffer
//

#if PY_BIG_ENDIAN
uint16_t reverse_16(uint16_t number) {
  return (number << 8) | (number >> 8);
}

uint32_t reverse_32(uint32_t number) {
  return ((number & 0xFF) << 24) | ((number & 0xFF00) << 8) | ((number & 0xFF0000) >> 8) | (number >> 24);
}
#endif

typedef struct {
    PyObject_HEAD
    char *ptr;  // Current read location in the buffer
    char *end;  // End of the buffer
    PyObject *source;  // The object that contains the buffer
} ReadBufferObject;

static PyTypeObject ReadBufferType;

static PyObject*
ReadBuffer_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    if (type != &ReadBufferType) {
        PyErr_SetString(PyExc_TypeError, "ReadBuffer should not be subclassed");
        return NULL;
    }

    ReadBufferObject *self = (ReadBufferObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->source = NULL;
        self->ptr = NULL;
        self->end = NULL;
    }
    return (PyObject *) self;
}

static int
ReadBuffer_init_internal(ReadBufferObject *self, PyObject *source) {
    if (!PyBytes_CheckExact(source)) {
        PyErr_SetString(PyExc_TypeError, "source must be a bytes object");
        return -1;
    }
    self->source = Py_NewRef(source);
    self->ptr = PyBytes_AS_STRING(source);
    self->end = self->ptr + PyBytes_GET_SIZE(source);
    return 0;
}

static PyObject*
ReadBuffer_internal(PyObject *source) {
    ReadBufferObject *self = (ReadBufferObject *)ReadBufferType.tp_alloc(&ReadBufferType, 0);
    if (self == NULL)
        return NULL;
    self->ptr = NULL;
    self->end = NULL;
    self->source = NULL;
    if (ReadBuffer_init_internal(self, source) == -1) {
        Py_DECREF(self);
        return NULL;
    }
    return (PyObject *)self;
}

static int
ReadBuffer_init(ReadBufferObject *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"source", NULL};
    PyObject *source = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &source))
        return -1;

    return ReadBuffer_init_internal(self, source);
}

static void
ReadBuffer_dealloc(ReadBufferObject *self)
{
    Py_CLEAR(self->source);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyMethodDef ReadBuffer_methods[] = {
    {NULL}  /* Sentinel */
};

static PyTypeObject ReadBufferType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "ReadBuffer",
    .tp_doc = PyDoc_STR("Mypy cache buffer objects"),
    .tp_basicsize = sizeof(ReadBufferObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = ReadBuffer_new,
    .tp_init = (initproc) ReadBuffer_init,
    .tp_dealloc = (destructor) ReadBuffer_dealloc,
    .tp_methods = ReadBuffer_methods,
};

//
// WriteBuffer
//

typedef struct {
    PyObject_HEAD
    char *buf;  // Beginning of the buffer
    char *ptr;  // Current write location in the buffer
    char *end;  // End of the buffer
} WriteBufferObject;

static PyTypeObject WriteBufferType;

static PyObject*
WriteBuffer_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    if (type != &WriteBufferType) {
        PyErr_SetString(PyExc_TypeError, "WriteBuffer cannot be subclassed");
        return NULL;
    }

    WriteBufferObject *self = (WriteBufferObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->buf = NULL;
        self->ptr = NULL;
        self->end = NULL;
    }
    return (PyObject *)self;
}

static int
WriteBuffer_init_internal(WriteBufferObject *self) {
    Py_ssize_t size = START_SIZE;
    self->buf = PyMem_Malloc(size + 1);
    if (self->buf == NULL) {
        PyErr_NoMemory();
        return -1;
    }
    self->ptr = self->buf;
    self->end = self->buf + size;
    return 0;
}

static PyObject*
WriteBuffer_internal(void) {
    WriteBufferObject *self = (WriteBufferObject *)WriteBufferType.tp_alloc(&WriteBufferType, 0);
    if (self == NULL)
        return NULL;
    self->buf = NULL;
    self->ptr = NULL;
    self->end = NULL;
    if (WriteBuffer_init_internal(self) == -1) {
        Py_DECREF(self);
        return NULL;
    }
    return (PyObject *)self;
}

static int
WriteBuffer_init(WriteBufferObject *self, PyObject *args, PyObject *kwds)
{
    if (!PyArg_ParseTuple(args, "")) {
        return -1;
    }

    if (kwds != NULL && PyDict_Size(kwds) > 0) {
        PyErr_SetString(PyExc_TypeError,
                        "WriteBuffer() takes no keyword arguments");
        return -1;
    }

    return WriteBuffer_init_internal(self);
}

static void
WriteBuffer_dealloc(WriteBufferObject *self)
{
    PyMem_Free(self->buf);
    self->buf = NULL;
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject*
WriteBuffer_getvalue_internal(PyObject *self)
{
    WriteBufferObject *obj = (WriteBufferObject *)self;
    return PyBytes_FromStringAndSize(obj->buf, obj->ptr - obj->buf);
}

static PyObject*
WriteBuffer_getvalue(WriteBufferObject *self, PyObject *Py_UNUSED(ignored))
{
    return PyBytes_FromStringAndSize(self->buf, self->ptr - self->buf);
}

static PyMethodDef WriteBuffer_methods[] = {
    {"getvalue", (PyCFunction) WriteBuffer_getvalue, METH_NOARGS,
     "Return the buffer content as bytes object"
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject WriteBufferType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "WriteBuffer",
    .tp_doc = PyDoc_STR("Mypy cache buffer objects"),
    .tp_basicsize = sizeof(WriteBufferObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = WriteBuffer_new,
    .tp_init = (initproc) WriteBuffer_init,
    .tp_dealloc = (destructor) WriteBuffer_dealloc,
    .tp_methods = WriteBuffer_methods,
};

// ----------

static inline char
_check_read_buffer(PyObject *data) {
    if (unlikely(Py_TYPE(data) != &ReadBufferType)) {
        PyErr_Format(
            PyExc_TypeError, "data must be a ReadBuffer object, got %s", Py_TYPE(data)->tp_name
        );
        return CPY_NONE_ERROR;
    }
    return CPY_NONE;
}

static inline char
_check_write_buffer(PyObject *data) {
    if (unlikely(Py_TYPE(data) != &WriteBufferType)) {
        PyErr_Format(
            PyExc_TypeError, "data must be a WriteBuffer object, got %s", Py_TYPE(data)->tp_name
        );
        return CPY_NONE_ERROR;
    }
    return CPY_NONE;
}

static inline char
_check_size(WriteBufferObject *data, Py_ssize_t need) {
    if (data->end - data->ptr >= need)
        return CPY_NONE;
    Py_ssize_t index = data->ptr - data->buf;
    Py_ssize_t target = index + need;
    Py_ssize_t size = data->end - data->buf;
    do {
        size *= 2;
    } while (target >= size);
    data->buf = PyMem_Realloc(data->buf, size);
    if (unlikely(data->buf == NULL)) {
        PyErr_NoMemory();
        return CPY_NONE_ERROR;
    }
    data->ptr = data->buf + index;
    data->end = data->buf + size;
    return CPY_NONE;
}

static inline char
_check_read(ReadBufferObject *data, Py_ssize_t need) {
    if (unlikely((data->end - data->ptr) < need)) {
        PyErr_SetString(PyExc_ValueError, "reading past the buffer end");
        return CPY_NONE_ERROR;
    }
    return CPY_NONE;
}

/*
bool format: single byte
    \x00 - False
    \x01 - True
*/

static char
read_bool_internal(PyObject *data) {
    _CHECK_READ(data, 1, CPY_BOOL_ERROR)
    char res;
    _READ(&res, data, char);
    if (unlikely((res != 0) & (res != 1))) {
        PyErr_SetString(PyExc_ValueError, "invalid bool value");
        return CPY_BOOL_ERROR;
    }
    return res;
}

static PyObject*
read_bool(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 1)) {
        PyErr_Format(PyExc_TypeError,
                     "read_bool() takes exactly 1 argument (%zu given)", nargs);
        return NULL;
    }
    PyObject *data = args[0];
    _CHECK_READ_BUFFER(data, NULL)
    char res = read_bool_internal(data);
    if (unlikely(res == CPY_BOOL_ERROR))
        return NULL;
    PyObject *retval = res ? Py_True : Py_False;
    Py_INCREF(retval);
    return retval;
}

static char
write_bool_internal(PyObject *data, char value) {
    _CHECK_WRITE(data, 1)
    _WRITE(data, char, value);
    return CPY_NONE;
}

static PyObject*
write_bool(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 2)) {
        PyErr_Format(PyExc_TypeError,
                     "write_bool() takes exactly 2 arguments (%zu given)", nargs);
        return NULL;
    }
    PyObject *data = args[0];
    PyObject *value = args[1];
    _CHECK_WRITE_BUFFER(data, NULL)
    if (unlikely(!PyBool_Check(value))) {
        PyErr_SetString(PyExc_TypeError, "value must be a bool");
        return NULL;
    }
    if (unlikely(write_bool_internal(data, Py_IsTrue(value)) == CPY_NONE_ERROR)) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

/*
str format: size as int (see below) followed by UTF-8 bytes
*/

static inline CPyTagged
_read_short_int(PyObject *data, uint8_t first) {
    uint8_t second;
    uint16_t two_more;
    if ((first & TWO_BYTES_INT_BIT) == 0) {
       // Note we use tagged ints since this function can return an error.
       return ((Py_ssize_t)(first >> 1) + MIN_ONE_BYTE_INT) << 1;
    }
    if ((first & FOUR_BYTES_INT_BIT) == 0) {
       _CHECK_READ(data, 1, CPY_INT_TAG)
       _READ(&second, data, uint8_t);
       return ((((Py_ssize_t)second) << 6) + (Py_ssize_t)(first >> 2) + MIN_TWO_BYTES_INT) << 1;
    }
    // The caller is responsible to verify this is called only for short ints.
    _CHECK_READ(data, 3, CPY_INT_TAG)
    // TODO: check if compilers emit optimal code for these two reads, and tweak if needed.
    _READ(&second, data, uint8_t);
    _READ(&two_more, data, uint16_t);
#if PY_BIG_ENDIAN
    two_more = reverse_16(two_more);
#endif
    Py_ssize_t higher = (((Py_ssize_t)two_more) << 13) + (((Py_ssize_t)second) << 5);
    return (higher + (Py_ssize_t)(first >> 3) + MIN_FOUR_BYTES_INT) << 1;
}

static PyObject*
read_str_internal(PyObject *data) {
    // Read string length.
    _CHECK_READ(data, 1, NULL)
    uint8_t first;
    _READ(&first, data, uint8_t);
    if (unlikely(first == LONG_INT_TRAILER)) {
        // Fail fast for invalid/tampered data.
        PyErr_SetString(PyExc_ValueError, "invalid str size");
        return NULL;
    }
    CPyTagged tagged_size = _read_short_int(data, first);
    if (tagged_size == CPY_INT_TAG)
        return NULL;
    if ((Py_ssize_t)tagged_size < 0) {
        // Fail fast for invalid/tampered data.
        PyErr_SetString(PyExc_ValueError, "invalid str size");
        return NULL;
    }
    Py_ssize_t size = tagged_size >> 1;
    // Read string content.
    char *ptr = ((ReadBufferObject *)data)->ptr;
    _CHECK_READ(data, size, NULL)
    PyObject *res = PyUnicode_FromStringAndSize(ptr, (Py_ssize_t)size);
    if (unlikely(res == NULL))
        return NULL;
    ((ReadBufferObject *)data)->ptr += size;
    return res;
}

static PyObject*
read_str(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 1)) {
        PyErr_Format(PyExc_TypeError,
                     "read_str() takes exactly 1 argument (%zu given)", nargs);
        return NULL;
    }
    PyObject *data = args[0];
    _CHECK_READ_BUFFER(data, NULL)
    return read_str_internal(data);
}

// The caller *must* check that real_value is within allowed range (29 bits).
static inline char
_write_short_int(PyObject *data, Py_ssize_t real_value) {
    if (real_value >= MIN_ONE_BYTE_INT && real_value <= MAX_ONE_BYTE_INT) {
        _CHECK_WRITE(data, 1)
        _WRITE(data, uint8_t, (uint8_t)(real_value - MIN_ONE_BYTE_INT) << 1);
    } else if (real_value >= MIN_TWO_BYTES_INT && real_value <= MAX_TWO_BYTES_INT) {
        _CHECK_WRITE(data, 2)
#if PY_BIG_ENDIAN
        uint16_t to_write = ((uint16_t)(real_value - MIN_TWO_BYTES_INT) << 2) | TWO_BYTES_INT_BIT;
        _WRITE(data, uint16_t, reverse_16(to_write));
#else
        _WRITE(data, uint16_t, ((uint16_t)(real_value - MIN_TWO_BYTES_INT) << 2) | TWO_BYTES_INT_BIT);
#endif
    } else {
        _CHECK_WRITE(data, 4)
#if PY_BIG_ENDIAN
        uint32_t to_write = ((uint32_t)(real_value - MIN_FOUR_BYTES_INT) << 3) | FOUR_BYTES_INT_TRAILER;
        _WRITE(data, uint32_t, reverse_32(to_write));
#else
        _WRITE(data, uint32_t, ((uint32_t)(real_value - MIN_FOUR_BYTES_INT) << 3) | FOUR_BYTES_INT_TRAILER);
#endif
    }
    return CPY_NONE;
}

static char
write_str_internal(PyObject *data, PyObject *value) {
    Py_ssize_t size;
    const char *chunk = PyUnicode_AsUTF8AndSize(value, &size);
    if (unlikely(chunk == NULL))
        return CPY_NONE_ERROR;

    // Write string length.
    if (likely(size >= MIN_FOUR_BYTES_INT && size <= MAX_FOUR_BYTES_INT)) {
        if (_write_short_int(data, size) == CPY_NONE_ERROR)
            return CPY_NONE_ERROR;
    } else {
        PyErr_SetString(PyExc_ValueError, "str too long to serialize");
        return CPY_NONE_ERROR;
    }
    // Write string content.
    _CHECK_WRITE(data, size)
    char *ptr = ((WriteBufferObject *)data)->ptr;
    memcpy(ptr, chunk, size);
    ((WriteBufferObject *)data)->ptr += size;
    return CPY_NONE;
}

static PyObject*
write_str(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 2)) {
        PyErr_Format(PyExc_TypeError,
                     "write_str() takes exactly 2 arguments (%zu given)", nargs);
        return NULL;
    }
    PyObject *data = args[0];
    PyObject *value = args[1];
    _CHECK_WRITE_BUFFER(data, NULL)
    if (unlikely(!PyUnicode_Check(value))) {
        PyErr_SetString(PyExc_TypeError, "value must be a str");
        return NULL;
    }
    if (unlikely(write_str_internal(data, value) == CPY_NONE_ERROR)) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

/*
bytes format: size as int (see below) followed by bytes
*/

static PyObject*
read_bytes_internal(PyObject *data) {
    // Read length.
    _CHECK_READ(data, 1, NULL)
    uint8_t first;
    _READ(&first, data, uint8_t);
    if (unlikely(first == LONG_INT_TRAILER)) {
        // Fail fast for invalid/tampered data.
        PyErr_SetString(PyExc_ValueError, "invalid bytes size");
        return NULL;
    }
    CPyTagged tagged_size = _read_short_int(data, first);
    if (tagged_size == CPY_INT_TAG)
        return NULL;
    if ((Py_ssize_t)tagged_size < 0) {
        // Fail fast for invalid/tampered data.
        PyErr_SetString(PyExc_ValueError, "invalid bytes size");
        return NULL;
    }
    Py_ssize_t size = tagged_size >> 1;
    // Read bytes content.
    char *ptr = ((ReadBufferObject *)data)->ptr;
    _CHECK_READ(data, size, NULL)
    PyObject *res = PyBytes_FromStringAndSize(ptr, (Py_ssize_t)size);
    if (unlikely(res == NULL))
        return NULL;
    ((ReadBufferObject *)data)->ptr += size;
    return res;
}

static PyObject*
read_bytes(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 1)) {
        PyErr_Format(PyExc_TypeError,
                     "read_bytes() takes exactly 1 argument (%zu given)", nargs);
        return NULL;
    }
    PyObject *data = args[0];
    _CHECK_READ_BUFFER(data, NULL)
    return read_bytes_internal(data);
}

static char
write_bytes_internal(PyObject *data, PyObject *value) {
    const char *chunk = PyBytes_AsString(value);
    if (unlikely(chunk == NULL))
        return CPY_NONE_ERROR;
    Py_ssize_t size = PyBytes_GET_SIZE(value);

    // Write length.
    if (likely(size >= MIN_FOUR_BYTES_INT && size <= MAX_FOUR_BYTES_INT)) {
        if (_write_short_int(data, size) == CPY_NONE_ERROR)
            return CPY_NONE_ERROR;
    } else {
        PyErr_SetString(PyExc_ValueError, "bytes too long to serialize");
        return CPY_NONE_ERROR;
    }
    // Write bytes content.
    _CHECK_WRITE(data, size)
    char *ptr = ((WriteBufferObject *)data)->ptr;
    memcpy(ptr, chunk, size);
    ((WriteBufferObject *)data)->ptr += size;
    return CPY_NONE;
}

static PyObject*
write_bytes(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 2)) {
        PyErr_Format(PyExc_TypeError,
                     "write_bytes() takes exactly 2 arguments (%zu given)", nargs);
        return NULL;
    }
    PyObject *data = args[0];
    PyObject *value = args[1];
    _CHECK_WRITE_BUFFER(data, NULL)
    if (unlikely(!PyBytes_Check(value))) {
        PyErr_SetString(PyExc_TypeError, "value must be a bytes object");
        return NULL;
    }
    if (unlikely(write_bytes_internal(data, value) == CPY_NONE_ERROR)) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

/*
float format:
    stored using PyFloat helpers in little-endian format.
*/

static double
read_float_internal(PyObject *data) {
    _CHECK_READ(data, 8, CPY_FLOAT_ERROR)
    char *ptr = ((ReadBufferObject *)data)->ptr;
    double res = PyFloat_Unpack8(ptr, 1);
    if (unlikely((res == -1.0) && PyErr_Occurred()))
        return CPY_FLOAT_ERROR;
    ((ReadBufferObject *)data)->ptr += 8;
    return res;
}

static PyObject*
read_float(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 1)) {
        PyErr_Format(PyExc_TypeError,
                     "read_float() takes exactly 1 argument (%zu given)", nargs);
        return NULL;
    }
    PyObject *data = args[0];
    _CHECK_READ_BUFFER(data, NULL)
    double retval = read_float_internal(data);
    if (unlikely(retval == CPY_FLOAT_ERROR && PyErr_Occurred())) {
        return NULL;
    }
    return PyFloat_FromDouble(retval);
}

static char
write_float_internal(PyObject *data, double value) {
    _CHECK_WRITE(data, 8)
    char *ptr = ((WriteBufferObject *)data)->ptr;
    int res = PyFloat_Pack8(value, ptr, 1);
    if (unlikely(res == -1))
        return CPY_NONE_ERROR;
    ((WriteBufferObject *)data)->ptr += 8;
    return CPY_NONE;
}

static PyObject*
write_float(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 2)) {
        PyErr_Format(PyExc_TypeError,
                     "write_float() takes exactly 2 arguments (%zu given)", nargs);
        return NULL;
    }
    PyObject *data = args[0];
    PyObject *value = args[1];
    _CHECK_WRITE_BUFFER(data, NULL)
    if (unlikely(!PyFloat_Check(value))) {
        PyErr_SetString(PyExc_TypeError, "value must be a float");
        return NULL;
    }
    if (unlikely(write_float_internal(data, PyFloat_AsDouble(value)) == CPY_NONE_ERROR)) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

/*
int format:
    one byte: last bit 0, 7 bits used
    two bytes: last two bits 01, 14 bits used
    four bytes: last three bits 011, 29 bits used
    everything else: 00001111 followed by serialized string representation

Note: for fixed size formats we skew ranges towards more positive values,
since negative integers are much more rare.
*/

static CPyTagged
read_int_internal(PyObject *data) {
    _CHECK_READ(data, 1, CPY_INT_TAG)

    uint8_t first;
    _READ(&first, data, uint8_t);
    if (likely(first != LONG_INT_TRAILER)) {
        return _read_short_int(data, first);
    }

    // Long integer encoding -- byte length and sign, followed by a byte array.

    // Read byte length and sign.
    _CHECK_READ(data, 1, CPY_INT_TAG)
    _READ(&first, data, uint8_t);
    Py_ssize_t size_and_sign = _read_short_int(data, first);
    if (size_and_sign == CPY_INT_TAG)
        return CPY_INT_TAG;
    if ((Py_ssize_t)size_and_sign < 0) {
        PyErr_SetString(PyExc_ValueError, "invalid int data");
        return CPY_INT_TAG;
    }
    bool sign = (size_and_sign >> 1) & 1;
    Py_ssize_t size = size_and_sign >> 2;

    // Construct an int object from the byte array.
    _CHECK_READ(data, size, CPY_INT_TAG)
    char *ptr = ((ReadBufferObject *)data)->ptr;
    PyObject *num = _PyLong_FromByteArray((unsigned char *)ptr, size, 1, 0);
    if (num == NULL)
        return CPY_INT_TAG;
    ((ReadBufferObject *)data)->ptr += size;
    if (sign) {
        PyObject *old = num;
        num = PyNumber_Negative(old);
        Py_DECREF(old);
        if (num == NULL) {
            return CPY_INT_TAG;
        }
    }
    return CPyTagged_StealFromObject(num);
}

static PyObject*
read_int(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 1)) {
        PyErr_Format(PyExc_TypeError,
                     "read_int() takes exactly 1 argument (%zu given)", nargs);
        return NULL;
    }
    PyObject *data = args[0];
    _CHECK_READ_BUFFER(data, NULL)
    CPyTagged retval = read_int_internal(data);
    if (unlikely(retval == CPY_INT_TAG)) {
        return NULL;
    }
    return CPyTagged_StealAsObject(retval);
}


static inline int hex_to_int(char c) {
    if (c >= '0' && c <= '9')
        return c - '0';
    else if (c >= 'a' && c <= 'f')
        return c - 'a' + 10;
    else
        return c - 'A' + 10;  // Assume valid hex digit
}

static inline char
_write_long_int(PyObject *data, CPyTagged value) {
    _CHECK_WRITE(data, 1)
    _WRITE(data, uint8_t, LONG_INT_TRAILER);

    PyObject *hex_str = NULL;
    PyObject* int_value = CPyTagged_AsObject(value);
    if (unlikely(int_value == NULL))
        goto error;

    hex_str = PyNumber_ToBase(int_value, 16);
    if (hex_str == NULL)
        goto error;
    Py_DECREF(int_value);
    int_value = NULL;

    const char *str = PyUnicode_AsUTF8(hex_str);
    if (str == NULL)
        goto error;
    Py_ssize_t len = strlen(str);
    bool neg;
    if (str[0] == '-') {
        str++;
        len--;
        neg = true;
    } else {
        neg = false;
    }
    // Skip the 0x hex prefix.
    str += 2;
    len -= 2;

    // Write bytes encoded length and sign.
    Py_ssize_t size = (len + 1) / 2;
    Py_ssize_t encoded_size = (size << 1) | neg;
    if (encoded_size <= MAX_FOUR_BYTES_INT) {
        if (_write_short_int(data, encoded_size) == CPY_NONE_ERROR)
            goto error;
    } else {
        PyErr_SetString(PyExc_ValueError, "int too long to serialize");
        goto error;
    }

    // Write absolute integer value as byte array in a variable-length little endian format.
    Py_ssize_t i;
    for (i = len; i > 1; i -= 2) {
        if (write_tag_internal(
                data, hex_to_int(str[i - 1]) | (hex_to_int(str[i - 2]) << 4)) == CPY_NONE_ERROR)
            goto error;
    }
    // The final byte may correspond to only one hex digit.
    if (i == 1) {
        if (write_tag_internal(data, hex_to_int(str[i - 1])) == CPY_NONE_ERROR)
            goto error;
    }

    Py_DECREF(hex_str);
    return CPY_NONE;

  error:

    Py_XDECREF(int_value);
    Py_XDECREF(hex_str);
    return CPY_NONE_ERROR;
}

static char
write_int_internal(PyObject *data, CPyTagged value) {
    if (likely((value & CPY_INT_TAG) == 0)) {
        Py_ssize_t real_value = CPyTagged_ShortAsSsize_t(value);
        if (likely(real_value >= MIN_FOUR_BYTES_INT && real_value <= MAX_FOUR_BYTES_INT)) {
            return _write_short_int(data, real_value);
        } else {
            return _write_long_int(data, value);
        }
    } else {
        return _write_long_int(data, value);
    }
}

static PyObject*
write_int(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 2)) {
        PyErr_Format(PyExc_TypeError,
                     "write_int() takes exactly 2 arguments (%zu given)", nargs);
        return NULL;
    }
    PyObject *data = args[0];
    PyObject *value = args[1];
    _CHECK_WRITE_BUFFER(data, NULL)
    if (unlikely(!PyLong_Check(value))) {
        PyErr_SetString(PyExc_TypeError, "value must be an int");
        return NULL;
    }
    CPyTagged tagged_value = CPyTagged_BorrowFromObject(value);
    if (unlikely(write_int_internal(data, tagged_value) == CPY_NONE_ERROR)) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

/*
integer tag format (0 <= t <= 255):
    stored as a uint8_t
*/

static uint8_t
read_tag_internal(PyObject *data) {
    _CHECK_READ(data, 1, CPY_LL_UINT_ERROR)
    uint8_t ret;
    _READ(&ret, data, uint8_t);
    return ret;
}

static PyObject*
read_tag(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 1)) {
        PyErr_Format(PyExc_TypeError,
                     "read_tag() takes exactly 1 argument (%zu given)", nargs);
        return NULL;
    }
    PyObject *data = args[0];
    _CHECK_READ_BUFFER(data, NULL)
    uint8_t retval = read_tag_internal(data);
    if (unlikely(retval == CPY_LL_UINT_ERROR && PyErr_Occurred())) {
        return NULL;
    }
    return PyLong_FromLong(retval);
}

static char
write_tag_internal(PyObject *data, uint8_t value) {
    _CHECK_WRITE(data, 1)
    _WRITE(data, uint8_t, value);
    return CPY_NONE;
}

static PyObject*
write_tag(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 2)) {
        PyErr_Format(PyExc_TypeError,
                     "write_tag() takes exactly 2 arguments (%zu given)", nargs);
        return NULL;
    }
    PyObject *data = args[0];
    PyObject *value = args[1];
    _CHECK_WRITE_BUFFER(data, NULL)
    uint8_t unboxed = CPyLong_AsUInt8(value);
    if (unlikely(unboxed == CPY_LL_UINT_ERROR && PyErr_Occurred())) {
        CPy_TypeError("u8", value);
        return NULL;
    }
    if (unlikely(write_tag_internal(data, unboxed) == CPY_NONE_ERROR)) {
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// All tags must be kept in sync with cache.py, nodes.py, and types.py.
// Primitive types.
#define LITERAL_FALSE 0
#define LITERAL_TRUE 1
#define LITERAL_NONE 2
#define LITERAL_INT 3
#define LITERAL_STR 4
#define LITERAL_BYTES 5
#define LITERAL_FLOAT 6
#define LITERAL_COMPLEX 7

// Supported builtin collections.
#define LIST_GEN 20
#define LIST_INT 21
#define LIST_STR 22
#define LIST_BYTES 23
#define TUPLE_GEN 24
#define DICT_STR_GEN 30

// This is the smallest custom class tag.
#define MYPY_FILE 50

// Instance class has special formats.
#define INSTANCE 80
#define INSTANCE_SIMPLE 81
#define INSTANCE_GENERIC 82
#define INSTANCE_STR 83
#define INSTANCE_FUNCTION 84
#define INSTANCE_INT 85
#define INSTANCE_BOOL 86
#define INSTANCE_OBJECT 87

#define RESERVED 254
#define END_TAG 255

// Forward declaration.
static char _skip_object(PyObject *data, uint8_t tag);

static inline char
_skip(PyObject *data, Py_ssize_t size) {
    // We are careful about error conditions, so all
    // _skip_xxx() functions can return an error value.
    _CHECK_READ(data, size, CPY_NONE_ERROR)
    ((ReadBufferObject *)data)->ptr += size;
    return CPY_NONE;
}

static inline char
_skip_short_int(PyObject *data, uint8_t first) {
    if ((first & TWO_BYTES_INT_BIT) == 0)
       return CPY_NONE;
    if ((first & FOUR_BYTES_INT_BIT) == 0)
        return _skip(data, 1);
    return _skip(data, 3);
}

static inline char
_skip_int(PyObject *data) {
    _CHECK_READ(data, 1, CPY_NONE_ERROR)

    uint8_t first;
    _READ(&first, data, uint8_t);
    if (likely(first != LONG_INT_TRAILER)) {
        return _skip_short_int(data, first);
    }

    _CHECK_READ(data, 1, CPY_NONE_ERROR)
    _READ(&first, data, uint8_t);
    Py_ssize_t size_and_sign = _read_short_int(data, first);
    if (size_and_sign == CPY_INT_TAG)
        return CPY_NONE_ERROR;
    if ((Py_ssize_t)size_and_sign < 0) {
        PyErr_SetString(PyExc_ValueError, "invalid int data");
        return CPY_NONE_ERROR;
    }
    Py_ssize_t size = size_and_sign >> 2;
    return _skip(data, size);
}

// This is essentially a wrapper around _read_short_int() that makes
// sure the result is valid.
static inline Py_ssize_t
_read_size(PyObject *data) {
    _CHECK_READ(data, 1, -1)
    uint8_t first;
    _READ(&first, data, uint8_t);
    // We actually allow serializing lists/dicts with over 4 billion items,
    // but we don't really need to, fail with ValueError just in case.
    if (unlikely(first == LONG_INT_TRAILER)) {
        PyErr_SetString(PyExc_ValueError, "unsupported size");
        return -1;
    }
    CPyTagged tagged_size = _read_short_int(data, first);
    if (tagged_size == CPY_INT_TAG)
        return -1;
    if ((Py_ssize_t)tagged_size < 0) {
        PyErr_SetString(PyExc_ValueError, "invalid size");
        return -1;
    }
    Py_ssize_t size = tagged_size >> 1;
    return size;
}

static inline char
_skip_str_bytes(PyObject *data) {
    Py_ssize_t size = _read_size(data);
    if (size < 0)
        return CPY_NONE_ERROR;
    return _skip(data, size);
}

// List/dict logic should be kept in sync with mypy/cache.py
static inline char
_skip_list_gen(PyObject *data) {
    Py_ssize_t size = _read_size(data);
    if (size < 0)
        return CPY_NONE_ERROR;
    Py_ssize_t i;
    for (i = 0; i < size; i++) {
        uint8_t tag = read_tag_internal(data);
        if (unlikely(tag == CPY_LL_UINT_ERROR && PyErr_Occurred())) {
            return CPY_NONE_ERROR;
        }
        if (unlikely(_skip_object(data, tag) == CPY_NONE_ERROR))
            return CPY_NONE_ERROR;
    }
    return CPY_NONE;
}

static inline char
_skip_list_int(PyObject *data) {
    Py_ssize_t size = _read_size(data);
    if (size < 0)
        return CPY_NONE_ERROR;
    Py_ssize_t i;
    for (i = 0; i < size; i++) {
        if (unlikely(_skip_int(data) == CPY_NONE_ERROR))
            return CPY_NONE_ERROR;
    }
    return CPY_NONE;
}

static inline char
_skip_list_str_bytes(PyObject *data) {
    Py_ssize_t size = _read_size(data);
    if (size < 0)
        return CPY_NONE_ERROR;
    Py_ssize_t i;
    for (i = 0; i < size; i++) {
        if (unlikely(_skip_str_bytes(data) == CPY_NONE_ERROR))
            return CPY_NONE_ERROR;
    }
    return CPY_NONE;
}

static inline char
_skip_dict_str_gen(PyObject *data) {
    Py_ssize_t size = _read_size(data);
    if (size < 0)
        return CPY_NONE_ERROR;
    Py_ssize_t i;
    for (i = 0; i < size; i++) {
        // Bare key followed by tagged value.
        if (unlikely(_skip_str_bytes(data) == CPY_NONE_ERROR))
            return CPY_NONE_ERROR;
        uint8_t tag = read_tag_internal(data);
        if (unlikely(tag == CPY_LL_UINT_ERROR && PyErr_Occurred())) {
            return CPY_NONE_ERROR;
        }
        if (unlikely(_skip_object(data, tag) == CPY_NONE_ERROR))
            return CPY_NONE_ERROR;
    }
    return CPY_NONE;
}

// Similar to mypy/cache.py, the convention is that the caller reads
// the opening tag for custom classes.
static inline char
_skip_class(PyObject *data) {
    while (1) {
        uint8_t tag = read_tag_internal(data);
        if (unlikely(tag == CPY_LL_UINT_ERROR && PyErr_Occurred())) {
            return CPY_NONE_ERROR;
        }
        if (tag == END_TAG) {
            return CPY_NONE;
        }
        if (unlikely(_skip_object(data, tag) == CPY_NONE_ERROR)) {
            return CPY_NONE_ERROR;
        }
    }
}

// Instance has special compact layout (as an important optimization).
static inline char
_skip_instance(PyObject *data) {
    uint8_t second_tag = read_tag_internal(data);
    if (unlikely(second_tag == CPY_LL_UINT_ERROR && PyErr_Occurred())) {
        return CPY_NONE_ERROR;
    }
    if (second_tag >= INSTANCE_STR && second_tag <= INSTANCE_OBJECT) {
        return CPY_NONE;
    }
    if (second_tag == INSTANCE_SIMPLE) {
        return _skip_str_bytes(data);
    }
    if (second_tag == INSTANCE_GENERIC) {
        return _skip_class(data);
    }
    PyErr_Format(PyExc_ValueError, "Unexpected instance tag: %d", second_tag);
    return CPY_NONE_ERROR;
}

// This is the main dispatch point. Branches are ordered manually
// based roughly on frequency in self-check.
static char
_skip_object(PyObject *data, uint8_t tag) {
    if (tag == LITERAL_STR || tag == LITERAL_BYTES)
        return _skip_str_bytes(data);
    if (tag == LITERAL_NONE || tag == LITERAL_FALSE || tag == LITERAL_TRUE)
        return CPY_NONE;
    if (tag == LIST_GEN || tag == TUPLE_GEN)
        return _skip_list_gen(data);
    if (tag == LITERAL_INT)
        return _skip_int(data);
    if (tag == INSTANCE)
        return _skip_instance(data);
    // We intentionally exclude MypyFile as a sanity check. Module symbols should
    // be always handled via cross_ref, and never appear in a symbol table as is.
    if (tag > MYPY_FILE && tag < RESERVED)
        return _skip_class(data);
    if (tag == LIST_INT)
        return _skip_list_int(data);
    if (tag == LIST_STR || tag == LIST_BYTES)
        return _skip_list_str_bytes(data);
    if (tag == DICT_STR_GEN)
        return _skip_dict_str_gen(data);
    if (tag == LITERAL_FLOAT)
        return _skip(data, 8);
    if (tag == LITERAL_COMPLEX)
        return _skip(data, 16);
    PyErr_Format(PyExc_ValueError, "Unsupported tag: %d", tag);
    return CPY_NONE_ERROR;
}

static PyObject*
extract_symbol_internal(PyObject *data) {
    char *ptr = ((ReadBufferObject *)data)->ptr;
    if (unlikely(_skip_class(data) == CPY_NONE_ERROR))
        return NULL;
    Py_ssize_t size = ((ReadBufferObject *)data)->ptr - ptr;
    PyObject *res = PyBytes_FromStringAndSize(ptr, size);
    if (unlikely(res == NULL))
        return NULL;
    return res;
}

static PyObject*
extract_symbol(PyObject *self, PyObject *const *args, size_t nargs) {
    if (unlikely(nargs != 1)) {
        PyErr_Format(PyExc_TypeError,
                     "extract_symbol() takes exactly 1 argument (%zu given)", nargs);
        return NULL;
    }
    PyObject *data = args[0];
    _CHECK_READ_BUFFER(data, NULL)
    return extract_symbol_internal(data);
}

static uint8_t
cache_version_internal(void) {
    return 0;
}

static PyObject*
cache_version(PyObject *self, PyObject *Py_UNUSED(ignored)) {
    return PyLong_FromLong(cache_version_internal());
}

static PyTypeObject *
ReadBuffer_type_internal(void) {
    return &ReadBufferType;  // Return borrowed reference
}

static PyTypeObject *
WriteBuffer_type_internal(void) {
    return &WriteBufferType;  // Return borrowed reference
};

static PyMethodDef librt_internal_module_methods[] = {
    {"write_bool", (PyCFunction)write_bool, METH_FASTCALL, PyDoc_STR("write a bool")},
    {"read_bool", (PyCFunction)read_bool, METH_FASTCALL, PyDoc_STR("read a bool")},
    {"write_str", (PyCFunction)write_str, METH_FASTCALL, PyDoc_STR("write a string")},
    {"read_str", (PyCFunction)read_str, METH_FASTCALL, PyDoc_STR("read a string")},
    {"write_bytes", (PyCFunction)write_bytes, METH_FASTCALL, PyDoc_STR("write bytes")},
    {"read_bytes", (PyCFunction)read_bytes, METH_FASTCALL, PyDoc_STR("read bytes")},
    {"write_float", (PyCFunction)write_float, METH_FASTCALL, PyDoc_STR("write a float")},
    {"read_float", (PyCFunction)read_float, METH_FASTCALL, PyDoc_STR("read a float")},
    {"write_int", (PyCFunction)write_int, METH_FASTCALL, PyDoc_STR("write an int")},
    {"read_int", (PyCFunction)read_int, METH_FASTCALL, PyDoc_STR("read an int")},
    {"write_tag", (PyCFunction)write_tag, METH_FASTCALL, PyDoc_STR("write a short int")},
    {"read_tag", (PyCFunction)read_tag, METH_FASTCALL, PyDoc_STR("read a short int")},
    {"cache_version", (PyCFunction)cache_version, METH_NOARGS, PyDoc_STR("cache format version")},
    {"extract_symbol", (PyCFunction)extract_symbol, METH_FASTCALL, PyDoc_STR("extract bytes for a mypy symbol")},
    {NULL, NULL, 0, NULL}
};

static int
NativeInternal_ABI_Version(void) {
    return LIBRT_INTERNAL_ABI_VERSION;
}

static int
NativeInternal_API_Version(void) {
    return LIBRT_INTERNAL_API_VERSION;
}

static int
librt_internal_module_exec(PyObject *m)
{
    if (PyType_Ready(&ReadBufferType) < 0) {
        return -1;
    }
    if (PyType_Ready(&WriteBufferType) < 0) {
        return -1;
    }
    if (PyModule_AddObjectRef(m, "ReadBuffer", (PyObject *) &ReadBufferType) < 0) {
        return -1;
    }
    if (PyModule_AddObjectRef(m, "WriteBuffer", (PyObject *) &WriteBufferType) < 0) {
        return -1;
    }

    // Export mypy internal C API, be careful with the order!
    static void *NativeInternal_API[LIBRT_INTERNAL_API_LEN] = {
        (void *)ReadBuffer_internal,
        (void *)WriteBuffer_internal,
        (void *)WriteBuffer_getvalue_internal,
        (void *)write_bool_internal,
        (void *)read_bool_internal,
        (void *)write_str_internal,
        (void *)read_str_internal,
        (void *)write_float_internal,
        (void *)read_float_internal,
        (void *)write_int_internal,
        (void *)read_int_internal,
        (void *)write_tag_internal,
        (void *)read_tag_internal,
        (void *)NativeInternal_ABI_Version,
        (void *)write_bytes_internal,
        (void *)read_bytes_internal,
        (void *)cache_version_internal,
        (void *)ReadBuffer_type_internal,
        (void *)WriteBuffer_type_internal,
        (void *)NativeInternal_API_Version,
        (void *)extract_symbol_internal
    };
    PyObject *c_api_object = PyCapsule_New((void *)NativeInternal_API, "librt.internal._C_API", NULL);
    if (PyModule_Add(m, "_C_API", c_api_object) < 0) {
        return -1;
    }
    return 0;
}

static PyModuleDef_Slot librt_internal_module_slots[] = {
    {Py_mod_exec, librt_internal_module_exec},
#ifdef Py_MOD_GIL_NOT_USED
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};

static PyModuleDef librt_internal_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "internal",
    .m_doc = "Mypy cache serialization utils",
    .m_size = 0,
    .m_methods = librt_internal_module_methods,
    .m_slots = librt_internal_module_slots,
};

PyMODINIT_FUNC
PyInit_internal(void)
{
    return PyModuleDef_Init(&librt_internal_module);
}
