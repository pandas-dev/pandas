#ifndef LIBRT_STRINGS_COMMON_H
#define LIBRT_STRINGS_COMMON_H

#include <Python.h>
#include <stdint.h>
#include <string.h>

// Byte-swap functions for endianness conversion (needed for both LE and BE operations)
#if defined(_MSC_VER)
#  include <stdlib.h>
#  define BSWAP16(x) _byteswap_ushort(x)
#  define BSWAP32(x) _byteswap_ulong(x)
#  define BSWAP64(x) _byteswap_uint64(x)
#elif defined(__GNUC__) || defined(__clang__)
#  define BSWAP16(x) __builtin_bswap16(x)
#  define BSWAP32(x) __builtin_bswap32(x)
#  define BSWAP64(x) __builtin_bswap64(x)
#else
// Fallback for other compilers (slower but portable)
static inline uint16_t BSWAP16(uint16_t x) {
    return (uint16_t)((x >> 8) | (x << 8));
}
static inline uint32_t BSWAP32(uint32_t x) {
    return ((x >> 24) & 0xFFU) |
           ((x >>  8) & 0xFF00U) |
           ((x <<  8) & 0xFF0000U) |
           ((x << 24) & 0xFF000000U);
}
static inline uint64_t BSWAP64(uint64_t x) {
    return ((x >> 56) & 0xFFULL) |
           ((x >> 40) & 0xFF00ULL) |
           ((x >> 24) & 0xFF0000ULL) |
           ((x >>  8) & 0xFF000000ULL) |
           ((x <<  8) & 0xFF00000000ULL) |
           ((x << 24) & 0xFF0000000000ULL) |
           ((x << 40) & 0xFF000000000000ULL) |
           ((x << 56) & 0xFF00000000000000ULL);
}
#endif

// Length of the default buffer embedded directly in a BytesWriter object
#define WRITER_EMBEDDED_BUF_LEN 256

typedef struct {
    PyObject_HEAD
    char *buf;  // Beginning of the buffer
    Py_ssize_t len;  // Current length (number of bytes written)
    Py_ssize_t capacity;  // Total capacity of the buffer
    char data[WRITER_EMBEDDED_BUF_LEN];  // Default buffer
} BytesWriterObject;

// Write a 16-bit signed integer in little-endian format to BytesWriter.
// NOTE: This does NOT check buffer capacity - caller must ensure space is available.
static inline void
BytesWriter_WriteI16LEUnsafe(BytesWriterObject *self, int16_t value) {
    // memcpy is reliably optimized to a single store by GCC, Clang, and MSVC
#if PY_BIG_ENDIAN
    uint16_t swapped = BSWAP16((uint16_t)value);
    memcpy(self->buf + self->len, &swapped, 2);
#else
    memcpy(self->buf + self->len, &value, 2);
#endif
    self->len += 2;
}

// Write a 16-bit signed integer in big-endian format to BytesWriter.
// NOTE: This does NOT check buffer capacity - caller must ensure space is available.
static inline void
BytesWriter_WriteI16BEUnsafe(BytesWriterObject *self, int16_t value) {
    // memcpy is reliably optimized to a single store by GCC, Clang, and MSVC
#if PY_BIG_ENDIAN
    memcpy(self->buf + self->len, &value, 2);
#else
    uint16_t swapped = BSWAP16((uint16_t)value);
    memcpy(self->buf + self->len, &swapped, 2);
#endif
    self->len += 2;
}

// Read a 16-bit signed integer in little-endian format from bytes.
// NOTE: This does NOT check bounds - caller must ensure valid index.
static inline int16_t
CPyBytes_ReadI16LEUnsafe(const unsigned char *data) {
    // memcpy is reliably optimized to a single load by GCC, Clang, and MSVC
    uint16_t value;
    memcpy(&value, data, 2);
#if PY_BIG_ENDIAN
    value = BSWAP16(value);
#endif
    return (int16_t)value;
}

// Read a 16-bit signed integer in big-endian format from bytes.
// NOTE: This does NOT check bounds - caller must ensure valid index.
static inline int16_t
CPyBytes_ReadI16BEUnsafe(const unsigned char *data) {
    // memcpy is reliably optimized to a single load by GCC, Clang, and MSVC
    uint16_t value;
    memcpy(&value, data, 2);
#if PY_BIG_ENDIAN
    // Already in big-endian format, no swap needed
#else
    value = BSWAP16(value);
#endif
    return (int16_t)value;
}

// Write a 32-bit signed integer in little-endian format to BytesWriter.
// NOTE: This does NOT check buffer capacity - caller must ensure space is available.
static inline void
BytesWriter_WriteI32LEUnsafe(BytesWriterObject *self, int32_t value) {
    // memcpy is reliably optimized to a single store by GCC, Clang, and MSVC
#if PY_BIG_ENDIAN
    uint32_t swapped = BSWAP32((uint32_t)value);
    memcpy(self->buf + self->len, &swapped, 4);
#else
    memcpy(self->buf + self->len, &value, 4);
#endif
    self->len += 4;
}

// Write a 32-bit signed integer in big-endian format to BytesWriter.
// NOTE: This does NOT check buffer capacity - caller must ensure space is available.
static inline void
BytesWriter_WriteI32BEUnsafe(BytesWriterObject *self, int32_t value) {
    // memcpy is reliably optimized to a single store by GCC, Clang, and MSVC
#if PY_BIG_ENDIAN
    memcpy(self->buf + self->len, &value, 4);
#else
    uint32_t swapped = BSWAP32((uint32_t)value);
    memcpy(self->buf + self->len, &swapped, 4);
#endif
    self->len += 4;
}

// Read a 32-bit signed integer in little-endian format from bytes.
// NOTE: This does NOT check bounds - caller must ensure valid index.
static inline int32_t
CPyBytes_ReadI32LEUnsafe(const unsigned char *data) {
    // memcpy is reliably optimized to a single load by GCC, Clang, and MSVC
    uint32_t value;
    memcpy(&value, data, 4);
#if PY_BIG_ENDIAN
    value = BSWAP32(value);
#endif
    return (int32_t)value;
}

// Read a 32-bit signed integer in big-endian format from bytes.
// NOTE: This does NOT check bounds - caller must ensure valid index.
static inline int32_t
CPyBytes_ReadI32BEUnsafe(const unsigned char *data) {
    // memcpy is reliably optimized to a single load by GCC, Clang, and MSVC
    uint32_t value;
    memcpy(&value, data, 4);
#if PY_BIG_ENDIAN
    // Already in big-endian format, no swap needed
#else
    value = BSWAP32(value);
#endif
    return (int32_t)value;
}

// Write a 64-bit signed integer in little-endian format to BytesWriter.
// NOTE: This does NOT check buffer capacity - caller must ensure space is available.
static inline void
BytesWriter_WriteI64LEUnsafe(BytesWriterObject *self, int64_t value) {
    // memcpy is reliably optimized to a single store by GCC, Clang, and MSVC
#if PY_BIG_ENDIAN
    uint64_t swapped = BSWAP64((uint64_t)value);
    memcpy(self->buf + self->len, &swapped, 8);
#else
    memcpy(self->buf + self->len, &value, 8);
#endif
    self->len += 8;
}

// Write a 64-bit signed integer in big-endian format to BytesWriter.
// NOTE: This does NOT check buffer capacity - caller must ensure space is available.
static inline void
BytesWriter_WriteI64BEUnsafe(BytesWriterObject *self, int64_t value) {
    // memcpy is reliably optimized to a single store by GCC, Clang, and MSVC
#if PY_BIG_ENDIAN
    memcpy(self->buf + self->len, &value, 8);
#else
    uint64_t swapped = BSWAP64((uint64_t)value);
    memcpy(self->buf + self->len, &swapped, 8);
#endif
    self->len += 8;
}

// Read a 64-bit signed integer in little-endian format from bytes.
// NOTE: This does NOT check bounds - caller must ensure valid index.
static inline int64_t
CPyBytes_ReadI64LEUnsafe(const unsigned char *data) {
    // memcpy is reliably optimized to a single load by GCC, Clang, and MSVC
    uint64_t value;
    memcpy(&value, data, 8);
#if PY_BIG_ENDIAN
    value = BSWAP64(value);
#endif
    return (int64_t)value;
}

// Read a 64-bit signed integer in big-endian format from bytes.
// NOTE: This does NOT check bounds - caller must ensure valid index.
static inline int64_t
CPyBytes_ReadI64BEUnsafe(const unsigned char *data) {
    // memcpy is reliably optimized to a single load by GCC, Clang, and MSVC
    uint64_t value;
    memcpy(&value, data, 8);
#if PY_BIG_ENDIAN
    // Already in big-endian format, no swap needed
#else
    value = BSWAP64(value);
#endif
    return (int64_t)value;
}

// Write a 32-bit float in little-endian format to BytesWriter.
// NOTE: This does NOT check buffer capacity - caller must ensure space is available.
static inline void
BytesWriter_WriteF32LEUnsafe(BytesWriterObject *self, float value) {
    // memcpy is reliably optimized to a single store by GCC, Clang, and MSVC
#if PY_BIG_ENDIAN
    uint32_t bits;
    memcpy(&bits, &value, 4);
    bits = BSWAP32(bits);
    memcpy(self->buf + self->len, &bits, 4);
#else
    memcpy(self->buf + self->len, &value, 4);
#endif
    self->len += 4;
}

// Write a 32-bit float in big-endian format to BytesWriter.
// NOTE: This does NOT check buffer capacity - caller must ensure space is available.
static inline void
BytesWriter_WriteF32BEUnsafe(BytesWriterObject *self, float value) {
    // memcpy is reliably optimized to a single store by GCC, Clang, and MSVC
#if PY_BIG_ENDIAN
    memcpy(self->buf + self->len, &value, 4);
#else
    uint32_t bits;
    memcpy(&bits, &value, 4);
    bits = BSWAP32(bits);
    memcpy(self->buf + self->len, &bits, 4);
#endif
    self->len += 4;
}

// Read a 32-bit float in little-endian format from bytes.
// NOTE: This does NOT check bounds - caller must ensure valid index.
static inline float
CPyBytes_ReadF32LEUnsafe(const unsigned char *data) {
    // memcpy is reliably optimized to a single load by GCC, Clang, and MSVC
    float value;
#if PY_BIG_ENDIAN
    uint32_t bits;
    memcpy(&bits, data, 4);
    bits = BSWAP32(bits);
    memcpy(&value, &bits, 4);
#else
    memcpy(&value, data, 4);
#endif
    return value;
}

// Read a 32-bit float in big-endian format from bytes.
// NOTE: This does NOT check bounds - caller must ensure valid index.
static inline float
CPyBytes_ReadF32BEUnsafe(const unsigned char *data) {
    // memcpy is reliably optimized to a single load by GCC, Clang, and MSVC
    float value;
#if PY_BIG_ENDIAN
    memcpy(&value, data, 4);
#else
    uint32_t bits;
    memcpy(&bits, data, 4);
    bits = BSWAP32(bits);
    memcpy(&value, &bits, 4);
#endif
    return value;
}

// Write a 64-bit float (double) in little-endian format to BytesWriter.
// NOTE: This does NOT check buffer capacity - caller must ensure space is available.
static inline void
BytesWriter_WriteF64LEUnsafe(BytesWriterObject *self, double value) {
    // memcpy is reliably optimized to a single store by GCC, Clang, and MSVC
#if PY_BIG_ENDIAN
    uint64_t bits;
    memcpy(&bits, &value, 8);
    bits = BSWAP64(bits);
    memcpy(self->buf + self->len, &bits, 8);
#else
    memcpy(self->buf + self->len, &value, 8);
#endif
    self->len += 8;
}

// Write a 64-bit float (double) in big-endian format to BytesWriter.
// NOTE: This does NOT check buffer capacity - caller must ensure space is available.
static inline void
BytesWriter_WriteF64BEUnsafe(BytesWriterObject *self, double value) {
    // memcpy is reliably optimized to a single store by GCC, Clang, and MSVC
#if PY_BIG_ENDIAN
    memcpy(self->buf + self->len, &value, 8);
#else
    uint64_t bits;
    memcpy(&bits, &value, 8);
    bits = BSWAP64(bits);
    memcpy(self->buf + self->len, &bits, 8);
#endif
    self->len += 8;
}

// Read a 64-bit float (double) in little-endian format from bytes.
// NOTE: This does NOT check bounds - caller must ensure valid index.
static inline double
CPyBytes_ReadF64LEUnsafe(const unsigned char *data) {
    // memcpy is reliably optimized to a single load by GCC, Clang, and MSVC
    double value;
#if PY_BIG_ENDIAN
    uint64_t bits;
    memcpy(&bits, data, 8);
    bits = BSWAP64(bits);
    memcpy(&value, &bits, 8);
#else
    memcpy(&value, data, 8);
#endif
    return value;
}

// Read a 64-bit float (double) in big-endian format from bytes.
// NOTE: This does NOT check bounds - caller must ensure valid index.
static inline double
CPyBytes_ReadF64BEUnsafe(const unsigned char *data) {
    // memcpy is reliably optimized to a single load by GCC, Clang, and MSVC
    double value;
#if PY_BIG_ENDIAN
    memcpy(&value, data, 8);
#else
    uint64_t bits;
    memcpy(&bits, data, 8);
    bits = BSWAP64(bits);
    memcpy(&value, &bits, 8);
#endif
    return value;
}

#endif  // LIBRT_STRINGS_COMMON_H
