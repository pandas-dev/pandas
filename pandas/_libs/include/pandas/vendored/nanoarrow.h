// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#ifndef NANOARROW_BUILD_ID_H_INCLUDED
#define NANOARROW_BUILD_ID_H_INCLUDED

#define NANOARROW_VERSION_MAJOR 0
#define NANOARROW_VERSION_MINOR 3
#define NANOARROW_VERSION_PATCH 0
#define NANOARROW_VERSION "0.3.0-SNAPSHOT"

#define NANOARROW_VERSION_INT                                        \
  (NANOARROW_VERSION_MAJOR * 10000 + NANOARROW_VERSION_MINOR * 100 + \
   NANOARROW_VERSION_PATCH)

// #define NANOARROW_NAMESPACE YourNamespaceHere

#endif
// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#ifndef NANOARROW_NANOARROW_TYPES_H_INCLUDED
#define NANOARROW_NANOARROW_TYPES_H_INCLUDED

#include <stdint.h>
#include <string.h>



#if defined(NANOARROW_DEBUG) && !defined(NANOARROW_PRINT_AND_DIE)
#include <stdio.h>
#include <stdlib.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Extra guard for versions of Arrow without the canonical guard
#ifndef ARROW_FLAG_DICTIONARY_ORDERED

/// \defgroup nanoarrow-arrow-cdata Arrow C Data interface
///
/// The Arrow C Data (https://arrow.apache.org/docs/format/CDataInterface.html)
/// and Arrow C Stream (https://arrow.apache.org/docs/format/CStreamInterface.html)
/// interfaces are part of the
/// Arrow Columnar Format specification
/// (https://arrow.apache.org/docs/format/Columnar.html). See the Arrow documentation for
/// documentation of these structures.
///
/// @{

#ifndef ARROW_C_DATA_INTERFACE
#define ARROW_C_DATA_INTERFACE

#define ARROW_FLAG_DICTIONARY_ORDERED 1
#define ARROW_FLAG_NULLABLE 2
#define ARROW_FLAG_MAP_KEYS_SORTED 4

struct ArrowSchema {
  // Array type description
  const char* format;
  const char* name;
  const char* metadata;
  int64_t flags;
  int64_t n_children;
  struct ArrowSchema** children;
  struct ArrowSchema* dictionary;

  // Release callback
  void (*release)(struct ArrowSchema*);
  // Opaque producer-specific data
  void* private_data;
};

struct ArrowArray {
  // Array data description
  int64_t length;
  int64_t null_count;
  int64_t offset;
  int64_t n_buffers;
  int64_t n_children;
  const void** buffers;
  struct ArrowArray** children;
  struct ArrowArray* dictionary;

  // Release callback
  void (*release)(struct ArrowArray*);
  // Opaque producer-specific data
  void* private_data;
};

#endif  // ARROW_C_DATA_INTERFACE

#ifndef ARROW_C_STREAM_INTERFACE
#define ARROW_C_STREAM_INTERFACE

struct ArrowArrayStream {
  // Callback to get the stream type
  // (will be the same for all arrays in the stream).
  //
  // Return value: 0 if successful, an `errno`-compatible error code otherwise.
  //
  // If successful, the ArrowSchema must be released independently from the stream.
  int (*get_schema)(struct ArrowArrayStream*, struct ArrowSchema* out);

  // Callback to get the next array
  // (if no error and the array is released, the stream has ended)
  //
  // Return value: 0 if successful, an `errno`-compatible error code otherwise.
  //
  // If successful, the ArrowArray must be released independently from the stream.
  int (*get_next)(struct ArrowArrayStream*, struct ArrowArray* out);

  // Callback to get optional detailed error information.
  // This must only be called if the last stream operation failed
  // with a non-0 return code.
  //
  // Return value: pointer to a null-terminated character array describing
  // the last error, or NULL if no description is available.
  //
  // The returned pointer is only valid until the next operation on this stream
  // (including release).
  const char* (*get_last_error)(struct ArrowArrayStream*);

  // Release callback: release the stream's own resources.
  // Note that arrays returned by `get_next` must be individually released.
  void (*release)(struct ArrowArrayStream*);

  // Opaque producer-specific data
  void* private_data;
};

#endif  // ARROW_C_STREAM_INTERFACE
#endif  // ARROW_FLAG_DICTIONARY_ORDERED

/// \brief Move the contents of src into dst and set src->release to NULL
static inline void ArrowSchemaMove(struct ArrowSchema* src, struct ArrowSchema* dst) {
  memcpy(dst, src, sizeof(struct ArrowSchema));
  src->release = NULL;
}

/// \brief Move the contents of src into dst and set src->release to NULL
static inline void ArrowArrayMove(struct ArrowArray* src, struct ArrowArray* dst) {
  memcpy(dst, src, sizeof(struct ArrowArray));
  src->release = NULL;
}

/// \brief Move the contents of src into dst and set src->release to NULL
static inline void ArrowArrayStreamMove(struct ArrowArrayStream* src,
                                        struct ArrowArrayStream* dst) {
  memcpy(dst, src, sizeof(struct ArrowArrayStream));
  src->release = NULL;
}

/// @}

// Utility macros
#define _NANOARROW_CONCAT(x, y) x##y
#define _NANOARROW_MAKE_NAME(x, y) _NANOARROW_CONCAT(x, y)

#define _NANOARROW_RETURN_NOT_OK_IMPL(NAME, EXPR) \
  do {                                            \
    const int NAME = (EXPR);                      \
    if (NAME) return NAME;                        \
  } while (0)

#define _NANOARROW_CHECK_RANGE(x_, min_, max_) \
  NANOARROW_RETURN_NOT_OK((x_ >= min_ && x_ <= max_) ? NANOARROW_OK : EINVAL)

#define _NANOARROW_CHECK_UPPER_LIMIT(x_, max_) \
  NANOARROW_RETURN_NOT_OK((x_ <= max_) ? NANOARROW_OK : EINVAL)

#if defined(NANOARROW_DEBUG)
#define _NANOARROW_RETURN_NOT_OK_WITH_ERROR_IMPL(NAME, EXPR, ERROR_PTR_EXPR, EXPR_STR) \
  do {                                                                                 \
    const int NAME = (EXPR);                                                           \
    if (NAME) {                                                                        \
      ArrowErrorSet((ERROR_PTR_EXPR), "%s failed with errno %d\n* %s:%d", EXPR_STR,    \
                    NAME, __FILE__, __LINE__);                                         \
      return NAME;                                                                     \
    }                                                                                  \
  } while (0)
#else
#define _NANOARROW_RETURN_NOT_OK_WITH_ERROR_IMPL(NAME, EXPR, ERROR_PTR_EXPR, EXPR_STR) \
  do {                                                                                 \
    const int NAME = (EXPR);                                                           \
    if (NAME) {                                                                        \
      ArrowErrorSet((ERROR_PTR_EXPR), "%s failed with errno %d", EXPR_STR, NAME);      \
      return NAME;                                                                     \
    }                                                                                  \
  } while (0)
#endif

/// \brief Return code for success.
/// \ingroup nanoarrow-errors
#define NANOARROW_OK 0

/// \brief Represents an errno-compatible error code
/// \ingroup nanoarrow-errors
typedef int ArrowErrorCode;

/// \brief Check the result of an expression and return it if not NANOARROW_OK
/// \ingroup nanoarrow-errors
#define NANOARROW_RETURN_NOT_OK(EXPR) \
  _NANOARROW_RETURN_NOT_OK_IMPL(_NANOARROW_MAKE_NAME(errno_status_, __COUNTER__), EXPR)

/// \brief Check the result of an expression and return it if not NANOARROW_OK,
/// adding an auto-generated message to an ArrowError.
/// \ingroup nanoarrow-errors
///
/// This macro is used to ensure that functions that accept an ArrowError
/// as input always set its message when returning an error code (e.g., when calling
/// a nanoarrow function that does *not* accept ArrowError).
#define NANOARROW_RETURN_NOT_OK_WITH_ERROR(EXPR, ERROR_EXPR) \
  _NANOARROW_RETURN_NOT_OK_WITH_ERROR_IMPL(                  \
      _NANOARROW_MAKE_NAME(errno_status_, __COUNTER__), EXPR, ERROR_EXPR, #EXPR)

#if defined(NANOARROW_DEBUG) && !defined(NANOARROW_PRINT_AND_DIE)
#define NANOARROW_PRINT_AND_DIE(VALUE, EXPR_STR)                                  \
  do {                                                                            \
    fprintf(stderr, "%s failed with errno %d\n* %s:%d\n", EXPR_STR, (int)(VALUE), \
            __FILE__, (int)__LINE__);                                             \
    abort();                                                                      \
  } while (0)
#endif

#if defined(NANOARROW_DEBUG)
#define _NANOARROW_ASSERT_OK_IMPL(NAME, EXPR, EXPR_STR) \
  do {                                                  \
    const int NAME = (EXPR);                            \
    if (NAME) NANOARROW_PRINT_AND_DIE(NAME, EXPR_STR);  \
  } while (0)

/// \brief Assert that an expression's value is NANOARROW_OK
/// \ingroup nanoarrow-errors
///
/// If nanoarrow was built in debug mode (i.e., defined(NANOARROW_DEBUG) is true),
/// print a message to stderr and abort. If nanoarrow was built in release mode,
/// this statement has no effect. You can customize fatal error behaviour
/// be defining the NANOARROW_PRINT_AND_DIE macro before including nanoarrow.h
/// This macro is provided as a convenience for users and is not used internally.
#define NANOARROW_ASSERT_OK(EXPR) \
  _NANOARROW_ASSERT_OK_IMPL(_NANOARROW_MAKE_NAME(errno_status_, __COUNTER__), EXPR, #EXPR)
#else
#define NANOARROW_ASSERT_OK(EXPR) EXPR
#endif

static char _ArrowIsLittleEndian(void) {
  uint32_t check = 1;
  char first_byte;
  memcpy(&first_byte, &check, sizeof(char));
  return first_byte;
}

/// \brief Arrow type enumerator
/// \ingroup nanoarrow-utils
///
/// These names are intended to map to the corresponding arrow::Type::type
/// enumerator; however, the numeric values are specifically not equal
/// (i.e., do not rely on numeric comparison).
enum ArrowType {
  NANOARROW_TYPE_UNINITIALIZED = 0,
  NANOARROW_TYPE_NA = 1,
  NANOARROW_TYPE_BOOL,
  NANOARROW_TYPE_UINT8,
  NANOARROW_TYPE_INT8,
  NANOARROW_TYPE_UINT16,
  NANOARROW_TYPE_INT16,
  NANOARROW_TYPE_UINT32,
  NANOARROW_TYPE_INT32,
  NANOARROW_TYPE_UINT64,
  NANOARROW_TYPE_INT64,
  NANOARROW_TYPE_HALF_FLOAT,
  NANOARROW_TYPE_FLOAT,
  NANOARROW_TYPE_DOUBLE,
  NANOARROW_TYPE_STRING,
  NANOARROW_TYPE_BINARY,
  NANOARROW_TYPE_FIXED_SIZE_BINARY,
  NANOARROW_TYPE_DATE32,
  NANOARROW_TYPE_DATE64,
  NANOARROW_TYPE_TIMESTAMP,
  NANOARROW_TYPE_TIME32,
  NANOARROW_TYPE_TIME64,
  NANOARROW_TYPE_INTERVAL_MONTHS,
  NANOARROW_TYPE_INTERVAL_DAY_TIME,
  NANOARROW_TYPE_DECIMAL128,
  NANOARROW_TYPE_DECIMAL256,
  NANOARROW_TYPE_LIST,
  NANOARROW_TYPE_STRUCT,
  NANOARROW_TYPE_SPARSE_UNION,
  NANOARROW_TYPE_DENSE_UNION,
  NANOARROW_TYPE_DICTIONARY,
  NANOARROW_TYPE_MAP,
  NANOARROW_TYPE_EXTENSION,
  NANOARROW_TYPE_FIXED_SIZE_LIST,
  NANOARROW_TYPE_DURATION,
  NANOARROW_TYPE_LARGE_STRING,
  NANOARROW_TYPE_LARGE_BINARY,
  NANOARROW_TYPE_LARGE_LIST,
  NANOARROW_TYPE_INTERVAL_MONTH_DAY_NANO
};

/// \brief Get a string value of an enum ArrowType value
/// \ingroup nanoarrow-utils
///
/// Returns NULL for invalid values for type
static inline const char* ArrowTypeString(enum ArrowType type);

static inline const char* ArrowTypeString(enum ArrowType type) {
  switch (type) {
    case NANOARROW_TYPE_NA:
      return "na";
    case NANOARROW_TYPE_BOOL:
      return "bool";
    case NANOARROW_TYPE_UINT8:
      return "uint8";
    case NANOARROW_TYPE_INT8:
      return "int8";
    case NANOARROW_TYPE_UINT16:
      return "uint16";
    case NANOARROW_TYPE_INT16:
      return "int16";
    case NANOARROW_TYPE_UINT32:
      return "uint32";
    case NANOARROW_TYPE_INT32:
      return "int32";
    case NANOARROW_TYPE_UINT64:
      return "uint64";
    case NANOARROW_TYPE_INT64:
      return "int64";
    case NANOARROW_TYPE_HALF_FLOAT:
      return "half_float";
    case NANOARROW_TYPE_FLOAT:
      return "float";
    case NANOARROW_TYPE_DOUBLE:
      return "double";
    case NANOARROW_TYPE_STRING:
      return "string";
    case NANOARROW_TYPE_BINARY:
      return "binary";
    case NANOARROW_TYPE_FIXED_SIZE_BINARY:
      return "fixed_size_binary";
    case NANOARROW_TYPE_DATE32:
      return "date32";
    case NANOARROW_TYPE_DATE64:
      return "date64";
    case NANOARROW_TYPE_TIMESTAMP:
      return "timestamp";
    case NANOARROW_TYPE_TIME32:
      return "time32";
    case NANOARROW_TYPE_TIME64:
      return "time64";
    case NANOARROW_TYPE_INTERVAL_MONTHS:
      return "interval_months";
    case NANOARROW_TYPE_INTERVAL_DAY_TIME:
      return "interval_day_time";
    case NANOARROW_TYPE_DECIMAL128:
      return "decimal128";
    case NANOARROW_TYPE_DECIMAL256:
      return "decimal256";
    case NANOARROW_TYPE_LIST:
      return "list";
    case NANOARROW_TYPE_STRUCT:
      return "struct";
    case NANOARROW_TYPE_SPARSE_UNION:
      return "sparse_union";
    case NANOARROW_TYPE_DENSE_UNION:
      return "dense_union";
    case NANOARROW_TYPE_DICTIONARY:
      return "dictionary";
    case NANOARROW_TYPE_MAP:
      return "map";
    case NANOARROW_TYPE_EXTENSION:
      return "extension";
    case NANOARROW_TYPE_FIXED_SIZE_LIST:
      return "fixed_size_list";
    case NANOARROW_TYPE_DURATION:
      return "duration";
    case NANOARROW_TYPE_LARGE_STRING:
      return "large_string";
    case NANOARROW_TYPE_LARGE_BINARY:
      return "large_binary";
    case NANOARROW_TYPE_LARGE_LIST:
      return "large_list";
    case NANOARROW_TYPE_INTERVAL_MONTH_DAY_NANO:
      return "interval_month_day_nano";
    default:
      return NULL;
  }
}

/// \brief Arrow time unit enumerator
/// \ingroup nanoarrow-utils
///
/// These names and values map to the corresponding arrow::TimeUnit::type
/// enumerator.
enum ArrowTimeUnit {
  NANOARROW_TIME_UNIT_SECOND = 0,
  NANOARROW_TIME_UNIT_MILLI = 1,
  NANOARROW_TIME_UNIT_MICRO = 2,
  NANOARROW_TIME_UNIT_NANO = 3
};

/// \brief Validation level enumerator
/// \ingroup nanoarrow-array
enum ArrowValidationLevel {
  /// \brief Do not validate buffer sizes or content.
  NANOARROW_VALIDATION_LEVEL_NONE = 0,

  /// \brief Validate buffer sizes that depend on array length but do not validate buffer
  /// sizes that depend on buffer data access.
  NANOARROW_VALIDATION_LEVEL_MINIMAL = 1,

  /// \brief Validate all buffer sizes, including those that require buffer data access,
  /// but do not perform any checks that are O(1) along the length of the buffers.
  NANOARROW_VALIDATION_LEVEL_DEFAULT = 2,

  /// \brief Validate all buffer sizes and all buffer content. This is useful in the
  /// context of untrusted input or input that may have been corrupted in transit.
  NANOARROW_VALIDATION_LEVEL_FULL = 3
};

/// \brief Get a string value of an enum ArrowTimeUnit value
/// \ingroup nanoarrow-utils
///
/// Returns NULL for invalid values for time_unit
static inline const char* ArrowTimeUnitString(enum ArrowTimeUnit time_unit);

static inline const char* ArrowTimeUnitString(enum ArrowTimeUnit time_unit) {
  switch (time_unit) {
    case NANOARROW_TIME_UNIT_SECOND:
      return "s";
    case NANOARROW_TIME_UNIT_MILLI:
      return "ms";
    case NANOARROW_TIME_UNIT_MICRO:
      return "us";
    case NANOARROW_TIME_UNIT_NANO:
      return "ns";
    default:
      return NULL;
  }
}

/// \brief Functional types of buffers as described in the Arrow Columnar Specification
/// \ingroup nanoarrow-array-view
enum ArrowBufferType {
  NANOARROW_BUFFER_TYPE_NONE,
  NANOARROW_BUFFER_TYPE_VALIDITY,
  NANOARROW_BUFFER_TYPE_TYPE_ID,
  NANOARROW_BUFFER_TYPE_UNION_OFFSET,
  NANOARROW_BUFFER_TYPE_DATA_OFFSET,
  NANOARROW_BUFFER_TYPE_DATA
};

/// \brief An non-owning view of a string
/// \ingroup nanoarrow-utils
struct ArrowStringView {
  /// \brief A pointer to the start of the string
  ///
  /// If size_bytes is 0, this value may be NULL.
  const char* data;

  /// \brief The size of the string in bytes,
  ///
  /// (Not including the null terminator.)
  int64_t size_bytes;
};

/// \brief Return a view of a const C string
/// \ingroup nanoarrow-utils
static inline struct ArrowStringView ArrowCharView(const char* value);

static inline struct ArrowStringView ArrowCharView(const char* value) {
  struct ArrowStringView out;

  out.data = value;
  if (value) {
    out.size_bytes = (int64_t)strlen(value);
  } else {
    out.size_bytes = 0;
  }

  return out;
}

union ArrowBufferViewData {
  const void* data;
  const int8_t* as_int8;
  const uint8_t* as_uint8;
  const int16_t* as_int16;
  const uint16_t* as_uint16;
  const int32_t* as_int32;
  const uint32_t* as_uint32;
  const int64_t* as_int64;
  const uint64_t* as_uint64;
  const double* as_double;
  const float* as_float;
  const char* as_char;
};

/// \brief An non-owning view of a buffer
/// \ingroup nanoarrow-utils
struct ArrowBufferView {
  /// \brief A pointer to the start of the buffer
  ///
  /// If size_bytes is 0, this value may be NULL.
  union ArrowBufferViewData data;

  /// \brief The size of the buffer in bytes
  int64_t size_bytes;
};

/// \brief Array buffer allocation and deallocation
/// \ingroup nanoarrow-buffer
///
/// Container for allocate, reallocate, and free methods that can be used
/// to customize allocation and deallocation of buffers when constructing
/// an ArrowArray.
struct ArrowBufferAllocator {
  /// \brief Reallocate a buffer or return NULL if it cannot be reallocated
  uint8_t* (*reallocate)(struct ArrowBufferAllocator* allocator, uint8_t* ptr,
                         int64_t old_size, int64_t new_size);

  /// \brief Deallocate a buffer allocated by this allocator
  void (*free)(struct ArrowBufferAllocator* allocator, uint8_t* ptr, int64_t size);

  /// \brief Opaque data specific to the allocator
  void* private_data;
};

/// \brief An owning mutable view of a buffer
/// \ingroup nanoarrow-buffer
struct ArrowBuffer {
  /// \brief A pointer to the start of the buffer
  ///
  /// If capacity_bytes is 0, this value may be NULL.
  uint8_t* data;

  /// \brief The size of the buffer in bytes
  int64_t size_bytes;

  /// \brief The capacity of the buffer in bytes
  int64_t capacity_bytes;

  /// \brief The allocator that will be used to reallocate and/or free the buffer
  struct ArrowBufferAllocator allocator;
};

/// \brief An owning mutable view of a bitmap
/// \ingroup nanoarrow-bitmap
struct ArrowBitmap {
  /// \brief An ArrowBuffer to hold the allocated memory
  struct ArrowBuffer buffer;

  /// \brief The number of bits that have been appended to the bitmap
  int64_t size_bits;
};

/// \brief A description of an arrangement of buffers
/// \ingroup nanoarrow-utils
///
/// Contains the minimum amount of information required to
/// calculate the size of each buffer in an ArrowArray knowing only
/// the length and offset of the array.
struct ArrowLayout {
  /// \brief The function of each buffer
  enum ArrowBufferType buffer_type[3];

  /// \brief The data type of each buffer
  enum ArrowType buffer_data_type[3];

  /// \brief The size of an element each buffer or 0 if this size is variable or unknown
  int64_t element_size_bits[3];

  /// \brief The number of elements in the child array per element in this array for a
  /// fixed-size list
  int64_t child_size_elements;
};

/// \brief A non-owning view of an ArrowArray
/// \ingroup nanoarrow-array-view
///
/// This data structure provides access to the values contained within
/// an ArrowArray with fields provided in a more readily-extractible
/// form. You can re-use an ArrowArrayView for multiple ArrowArrays
/// with the same storage type, use it to represent a hypothetical
/// ArrowArray that does not exist yet, or use it to validate the buffers
/// of a future ArrowArray.
struct ArrowArrayView {
  /// \brief The underlying ArrowArray or NULL if it has not been set or
  /// if the buffers in this ArrowArrayView are not backed by an ArrowArray.
  struct ArrowArray* array;

  /// \brief The number of elements from the physical start of the buffers.
  int64_t offset;

  /// \brief The number of elements in this view.
  int64_t length;

  /// \brief A cached null count or -1 to indicate that this value is unknown.
  int64_t null_count;

  /// \brief The type used to store values in this array
  ///
  /// This type represents only the minimum required information to
  /// extract values from the array buffers (e.g., for a Date32 array,
  /// this value will be NANOARROW_TYPE_INT32). For dictionary-encoded
  /// arrays, this will be the index type.
  enum ArrowType storage_type;

  /// \brief The buffer types, strides, and sizes of this Array's buffers
  struct ArrowLayout layout;

  /// \brief This Array's buffers as ArrowBufferView objects
  struct ArrowBufferView buffer_views[3];

  /// \brief The number of children of this view
  int64_t n_children;

  /// \brief Pointers to views of this array's children
  struct ArrowArrayView** children;

  /// \brief Pointer to a view of this array's dictionary
  struct ArrowArrayView* dictionary;

  /// \brief Union type id to child index mapping
  ///
  /// If storage_type is a union type, a 256-byte ArrowMalloc()ed buffer
  /// such that child_index == union_type_id_map[type_id] and
  /// type_id == union_type_id_map[128 + child_index]. This value may be
  /// NULL in the case where child_id == type_id.
  int8_t* union_type_id_map;
};

// Used as the private data member for ArrowArrays allocated here and accessed
// internally within inline ArrowArray* helpers.
struct ArrowArrayPrivateData {
  // Holder for the validity buffer (or first buffer for union types, which are
  // the only type whose first buffer is not a valdiity buffer)
  struct ArrowBitmap bitmap;

  // Holder for additional buffers as required
  struct ArrowBuffer buffers[2];

  // The array of pointers to buffers. This must be updated after a sequence
  // of appends to synchronize its values with the actual buffer addresses
  // (which may have ben reallocated uring that time)
  const void* buffer_data[3];

  // The storage data type, or NANOARROW_TYPE_UNINITIALIZED if unknown
  enum ArrowType storage_type;

  // The buffer arrangement for the storage type
  struct ArrowLayout layout;

  // Flag to indicate if there are non-sequence union type ids.
  // In the future this could be replaced with a type id<->child mapping
  // to support constructing unions in append mode where type_id != child_index
  int8_t union_type_id_is_child_index;
};

/// \brief A representation of an interval.
/// \ingroup nanoarrow-utils
struct ArrowInterval {
  /// \brief The type of interval being used
  enum ArrowType type;
  /// \brief The number of months represented by the interval
  int32_t months;
  /// \brief The number of days represented by the interval
  int32_t days;
  /// \brief The number of ms represented by the interval
  int32_t ms;
  /// \brief The number of ns represented by the interval
  int64_t ns;
};

/// \brief Zero initialize an Interval with a given unit
/// \ingroup nanoarrow-utils
static inline void ArrowIntervalInit(struct ArrowInterval* interval,
                                     enum ArrowType type) {
  memset(interval, 0, sizeof(struct ArrowInterval));
  interval->type = type;
}

/// \brief A representation of a fixed-precision decimal number
/// \ingroup nanoarrow-utils
///
/// This structure should be initialized with ArrowDecimalInit() once and
/// values set using ArrowDecimalSetInt(), ArrowDecimalSetBytes128(),
/// or ArrowDecimalSetBytes256().
struct ArrowDecimal {
  /// \brief An array of 64-bit integers of n_words length defined in native-endian order
  uint64_t words[4];

  /// \brief The number of significant digits this decimal number can represent
  int32_t precision;

  /// \brief The number of digits after the decimal point. This can be negative.
  int32_t scale;

  /// \brief The number of words in the words array
  int n_words;

  /// \brief Cached value used by the implementation
  int high_word_index;

  /// \brief Cached value used by the implementation
  int low_word_index;
};

/// \brief Initialize a decimal with a given set of type parameters
/// \ingroup nanoarrow-utils
static inline void ArrowDecimalInit(struct ArrowDecimal* decimal, int32_t bitwidth,
                                    int32_t precision, int32_t scale) {
  memset(decimal->words, 0, sizeof(decimal->words));
  decimal->precision = precision;
  decimal->scale = scale;
  decimal->n_words = bitwidth / 8 / sizeof(uint64_t);

  if (_ArrowIsLittleEndian()) {
    decimal->low_word_index = 0;
    decimal->high_word_index = decimal->n_words - 1;
  } else {
    decimal->low_word_index = decimal->n_words - 1;
    decimal->high_word_index = 0;
  }
}

/// \brief Get a signed integer value of a sufficiently small ArrowDecimal
///
/// This does not check if the decimal's precision sufficiently small to fit
/// within the signed 64-bit integer range (A precision less than or equal
/// to 18 is sufficiently small).
static inline int64_t ArrowDecimalGetIntUnsafe(struct ArrowDecimal* decimal) {
  return (int64_t)decimal->words[decimal->low_word_index];
}

/// \brief Copy the bytes of this decimal into a sufficiently large buffer
/// \ingroup nanoarrow-utils
static inline void ArrowDecimalGetBytes(struct ArrowDecimal* decimal, uint8_t* out) {
  memcpy(out, decimal->words, decimal->n_words * sizeof(uint64_t));
}

/// \brief Returns 1 if the value represented by decimal is >= 0 or -1 otherwise
/// \ingroup nanoarrow-utils
static inline int64_t ArrowDecimalSign(struct ArrowDecimal* decimal) {
  return 1 | ((int64_t)(decimal->words[decimal->high_word_index]) >> 63);
}

/// \brief Sets the integer value of this decimal
/// \ingroup nanoarrow-utils
static inline void ArrowDecimalSetInt(struct ArrowDecimal* decimal, int64_t value) {
  if (value < 0) {
    memset(decimal->words, 0xff, decimal->n_words * sizeof(uint64_t));
  } else {
    memset(decimal->words, 0, decimal->n_words * sizeof(uint64_t));
  }

  decimal->words[decimal->low_word_index] = value;
}

/// \brief Copy bytes from a buffer into this decimal
/// \ingroup nanoarrow-utils
static inline void ArrowDecimalSetBytes(struct ArrowDecimal* decimal,
                                        const uint8_t* value) {
  memcpy(decimal->words, value, decimal->n_words * sizeof(uint64_t));
}

#ifdef __cplusplus
}
#endif

#endif
// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#ifndef NANOARROW_H_INCLUDED
#define NANOARROW_H_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>



// If using CMake, optionally pass -DNANOARROW_NAMESPACE=MyNamespace which will set this
// define in nanoarrow_config.h. If not, you can optionally #define NANOARROW_NAMESPACE
// MyNamespace here.

// This section remaps the non-prefixed symbols to the prefixed symbols so that
// code written against this build can be used independent of the value of
// NANOARROW_NAMESPACE.
#ifdef NANOARROW_NAMESPACE
#define NANOARROW_CAT(A, B) A##B
#define NANOARROW_SYMBOL(A, B) NANOARROW_CAT(A, B)

#define ArrowNanoarrowVersion NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowNanoarrowVersion)
#define ArrowNanoarrowVersionInt \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowNanoarrowVersionInt)
#define ArrowErrorMessage NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowErrorMessage)
#define ArrowMalloc NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowMalloc)
#define ArrowRealloc NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowRealloc)
#define ArrowFree NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowFree)
#define ArrowBufferAllocatorDefault \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowBufferAllocatorDefault)
#define ArrowBufferDeallocator \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowBufferDeallocator)
#define ArrowErrorSet NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowErrorSet)
#define ArrowLayoutInit NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowLayoutInit)
#define ArrowSchemaInit NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaInit)
#define ArrowSchemaInitFromType \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaInitFromType)
#define ArrowSchemaSetType NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaSetType)
#define ArrowSchemaSetTypeStruct \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaSetTypeStruct)
#define ArrowSchemaSetTypeFixedSize \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaSetTypeFixedSize)
#define ArrowSchemaSetTypeDecimal \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaSetTypeDecimal)
#define ArrowSchemaSetTypeDateTime \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaSetTypeDateTime)
#define ArrowSchemaSetTypeUnion \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaSetTypeUnion)
#define ArrowSchemaDeepCopy NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaDeepCopy)
#define ArrowSchemaSetFormat NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaSetFormat)
#define ArrowSchemaSetName NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaSetName)
#define ArrowSchemaSetMetadata \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaSetMetadata)
#define ArrowSchemaAllocateChildren \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaAllocateChildren)
#define ArrowSchemaAllocateDictionary \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaAllocateDictionary)
#define ArrowMetadataReaderInit \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowMetadataReaderInit)
#define ArrowMetadataReaderRead \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowMetadataReaderRead)
#define ArrowMetadataSizeOf NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowMetadataSizeOf)
#define ArrowMetadataHasKey NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowMetadataHasKey)
#define ArrowMetadataGetValue NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowMetadataGetValue)
#define ArrowMetadataBuilderInit \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowMetadataBuilderInit)
#define ArrowMetadataBuilderAppend \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowMetadataBuilderAppend)
#define ArrowMetadataBuilderSet \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowMetadataBuilderSet)
#define ArrowMetadataBuilderRemove \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowMetadataBuilderRemove)
#define ArrowSchemaViewInit NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaViewInit)
#define ArrowSchemaToString NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowSchemaToString)
#define ArrowArrayInitFromType \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayInitFromType)
#define ArrowArrayInitFromSchema \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayInitFromSchema)
#define ArrowArrayInitFromArrayView \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayInitFromArrayView)
#define ArrowArrayInitFromArrayView \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayInitFromArrayView)
#define ArrowArrayAllocateChildren \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayAllocateChildren)
#define ArrowArrayAllocateDictionary \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayAllocateDictionary)
#define ArrowArraySetValidityBitmap \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArraySetValidityBitmap)
#define ArrowArraySetBuffer NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArraySetBuffer)
#define ArrowArrayReserve NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayReserve)
#define ArrowArrayFinishBuilding \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayFinishBuilding)
#define ArrowArrayFinishBuildingDefault \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayFinishBuildingDefault)
#define ArrowArrayViewInitFromType \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayViewInitFromType)
#define ArrowArrayViewInitFromSchema \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayViewInitFromSchema)
#define ArrowArrayViewAllocateChildren \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayViewAllocateChildren)
#define ArrowArrayViewAllocateDictionary \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayViewAllocateDictionary)
#define ArrowArrayViewSetLength \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayViewSetLength)
#define ArrowArrayViewSetArray \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayViewSetArray)
#define ArrowArrayViewSetArrayMinimal \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayViewSetArrayMinimal)
#define ArrowArrayViewValidate \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayViewValidate)
#define ArrowArrayViewReset NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowArrayViewReset)
#define ArrowBasicArrayStreamInit \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowBasicArrayStreamInit)
#define ArrowBasicArrayStreamSetArray \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowBasicArrayStreamSetArray)
#define ArrowBasicArrayStreamValidate \
  NANOARROW_SYMBOL(NANOARROW_NAMESPACE, ArrowBasicArrayStreamValidate)

#endif

#ifdef __cplusplus
extern "C" {
#endif

/// \defgroup nanoarrow Nanoarrow C library
///
/// Except where noted, objects are not thread-safe and clients should
/// take care to serialize accesses to methods.
///
/// Because this library is intended to be vendored, it provides full type
/// definitions and encourages clients to stack or statically allocate
/// where convenient.

/// \defgroup nanoarrow-malloc Memory management
///
/// Non-buffer members of a struct ArrowSchema and struct ArrowArray
/// must be allocated using ArrowMalloc() or ArrowRealloc() and freed
/// using ArrowFree() for schemas and arrays allocated here. Buffer members
/// are allocated using an ArrowBufferAllocator.
///
/// @{

/// \brief Allocate like malloc()
void* ArrowMalloc(int64_t size);

/// \brief Reallocate like realloc()
void* ArrowRealloc(void* ptr, int64_t size);

/// \brief Free a pointer allocated using ArrowMalloc() or ArrowRealloc().
void ArrowFree(void* ptr);

/// \brief Return the default allocator
///
/// The default allocator uses ArrowMalloc(), ArrowRealloc(), and
/// ArrowFree().
struct ArrowBufferAllocator ArrowBufferAllocatorDefault(void);

/// \brief Create a custom deallocator
///
/// Creates a buffer allocator with only a free method that can be used to
/// attach a custom deallocator to an ArrowBuffer. This may be used to
/// avoid copying an existing buffer that was not allocated using the
/// infrastructure provided here (e.g., by an R or Python object).
struct ArrowBufferAllocator ArrowBufferDeallocator(
    void (*custom_free)(struct ArrowBufferAllocator* allocator, uint8_t* ptr,
                        int64_t size),
    void* private_data);

/// @}

/// \defgroup nanoarrow-errors Error handling
///
/// Functions generally return an errno-compatible error code; functions that
/// need to communicate more verbose error information accept a pointer
/// to an ArrowError. This can be stack or statically allocated. The
/// content of the message is undefined unless an error code has been
/// returned. If a nanoarrow function is passed a non-null ArrowError pointer, the
/// ArrowError pointed to by the argument will be propagated with a
/// null-terminated error message. It is safe to pass a NULL ArrowError anywhere
/// in the nanoarrow API.
///
/// Except where documented, it is generally not safe to continue after a
/// function has returned a non-zero ArrowErrorCode. The NANOARROW_RETURN_NOT_OK and
/// NANOARROW_ASSERT_OK macros are provided to help propagate errors. C++ clients can use
/// the helpers provided in the nanoarrow.hpp header to facilitate using C++ idioms
/// for memory management and error propgagtion.
///
/// @{

/// \brief Error type containing a UTF-8 encoded message.
struct ArrowError {
  /// \brief A character buffer with space for an error message.
  char message[1024];
};

/// \brief Ensure an ArrowError is null-terminated by zeroing the first character.
///
/// If error is NULL, this function does nothing.
static inline void ArrowErrorInit(struct ArrowError* error) {
  if (error) {
    error->message[0] = '\0';
  }
}

/// \brief Set the contents of an error using printf syntax.
///
/// If error is NULL, this function does nothing and returns NANOARROW_OK.
ArrowErrorCode ArrowErrorSet(struct ArrowError* error, const char* fmt, ...);

/// \brief Get the contents of an error
///
/// If error is NULL, returns "", or returns the contents of the error message
/// otherwise.
const char* ArrowErrorMessage(struct ArrowError* error);

/// @}

/// \defgroup nanoarrow-utils Utility data structures
///
/// @{

/// \brief Return a version string in the form "major.minor.patch"
const char* ArrowNanoarrowVersion(void);

/// \brief Return an integer that can be used to compare versions sequentially
int ArrowNanoarrowVersionInt(void);

/// \brief Initialize a description of buffer arrangements from a storage type
void ArrowLayoutInit(struct ArrowLayout* layout, enum ArrowType storage_type);

/// \brief Create a string view from a null-terminated string
static inline struct ArrowStringView ArrowCharView(const char* value);

/// @}

/// \defgroup nanoarrow-schema Creating schemas
///
/// These functions allocate, copy, and destroy ArrowSchema structures
///
/// @{

/// \brief Initialize an ArrowSchema
///
/// Initializes the fields and release callback of schema_out. Caller
/// is responsible for calling the schema->release callback if
/// NANOARROW_OK is returned.
void ArrowSchemaInit(struct ArrowSchema* schema);

/// \brief Initialize an ArrowSchema from an ArrowType
///
/// A convenience constructor for that calls ArrowSchemaInit() and
/// ArrowSchemaSetType() for the common case of constructing an
/// unparameterized type. The caller is responsible for calling the schema->release
/// callback if NANOARROW_OK is returned.
ArrowErrorCode ArrowSchemaInitFromType(struct ArrowSchema* schema, enum ArrowType type);

/// \brief Get a human-readable summary of a Schema
///
/// Writes a summary of an ArrowSchema to out (up to n - 1 characters)
/// and returns the number of characters required for the output if
/// n were sufficiently large. If recursive is non-zero, the result will
/// also include children.
int64_t ArrowSchemaToString(struct ArrowSchema* schema, char* out, int64_t n,
                            char recursive);

/// \brief Set the format field of a schema from an ArrowType
///
/// Initializes the fields and release callback of schema_out. For
/// NANOARROW_TYPE_LIST, NANOARROW_TYPE_LARGE_LIST, and
/// NANOARROW_TYPE_MAP, the appropriate number of children are
/// allocated, initialized, and named; however, the caller must
/// ArrowSchemaSetType() on the preinitialized children. Schema must have been initialized
/// using ArrowSchemaInit() or ArrowSchemaDeepCopy().
ArrowErrorCode ArrowSchemaSetType(struct ArrowSchema* schema, enum ArrowType type);

/// \brief Set the format field and initialize children of a struct schema
///
/// The specified number of children are initialized; however, the caller is responsible
/// for calling ArrowSchemaSetType() and ArrowSchemaSetName() on each child.
/// Schema must have been initialized using ArrowSchemaInit() or ArrowSchemaDeepCopy().
ArrowErrorCode ArrowSchemaSetTypeStruct(struct ArrowSchema* schema, int64_t n_children);

/// \brief Set the format field of a fixed-size schema
///
/// Returns EINVAL for fixed_size <= 0 or for type that is not
/// NANOARROW_TYPE_FIXED_SIZE_BINARY or NANOARROW_TYPE_FIXED_SIZE_LIST.
/// For NANOARROW_TYPE_FIXED_SIZE_LIST, the appropriate number of children are
/// allocated, initialized, and named; however, the caller must
/// ArrowSchemaSetType() the first child. Schema must have been initialized using
/// ArrowSchemaInit() or ArrowSchemaDeepCopy().
ArrowErrorCode ArrowSchemaSetTypeFixedSize(struct ArrowSchema* schema,
                                           enum ArrowType type, int32_t fixed_size);

/// \brief Set the format field of a decimal schema
///
/// Returns EINVAL for scale <= 0 or for type that is not
/// NANOARROW_TYPE_DECIMAL128 or NANOARROW_TYPE_DECIMAL256. Schema must have been
/// initialized using ArrowSchemaInit() or ArrowSchemaDeepCopy().
ArrowErrorCode ArrowSchemaSetTypeDecimal(struct ArrowSchema* schema, enum ArrowType type,
                                         int32_t decimal_precision,
                                         int32_t decimal_scale);

/// \brief Set the format field of a time, timestamp, or duration schema
///
/// Returns EINVAL for type that is not
/// NANOARROW_TYPE_TIME32, NANOARROW_TYPE_TIME64,
/// NANOARROW_TYPE_TIMESTAMP, or NANOARROW_TYPE_DURATION. The
/// timezone parameter must be NULL for a non-timestamp type. Schema must have been
/// initialized using ArrowSchemaInit() or ArrowSchemaDeepCopy().
ArrowErrorCode ArrowSchemaSetTypeDateTime(struct ArrowSchema* schema, enum ArrowType type,
                                          enum ArrowTimeUnit time_unit,
                                          const char* timezone);

/// \brief Seet the format field of a union schema
///
/// Returns EINVAL for a type that is not NANOARROW_TYPE_DENSE_UNION
/// or NANOARROW_TYPE_SPARSE_UNION. The specified number of children are
/// allocated, and initialized.
ArrowErrorCode ArrowSchemaSetTypeUnion(struct ArrowSchema* schema, enum ArrowType type,
                                       int64_t n_children);

/// \brief Make a (recursive) copy of a schema
///
/// Allocates and copies fields of schema into schema_out.
ArrowErrorCode ArrowSchemaDeepCopy(struct ArrowSchema* schema,
                                   struct ArrowSchema* schema_out);

/// \brief Copy format into schema->format
///
/// schema must have been allocated using ArrowSchemaInitFromType() or
/// ArrowSchemaDeepCopy().
ArrowErrorCode ArrowSchemaSetFormat(struct ArrowSchema* schema, const char* format);

/// \brief Copy name into schema->name
///
/// schema must have been allocated using ArrowSchemaInitFromType() or
/// ArrowSchemaDeepCopy().
ArrowErrorCode ArrowSchemaSetName(struct ArrowSchema* schema, const char* name);

/// \brief Copy metadata into schema->metadata
///
/// schema must have been allocated using ArrowSchemaInitFromType() or
/// ArrowSchemaDeepCopy.
ArrowErrorCode ArrowSchemaSetMetadata(struct ArrowSchema* schema, const char* metadata);

/// \brief Allocate the schema->children array
///
/// Includes the memory for each child struct ArrowSchema.
/// schema must have been allocated using ArrowSchemaInitFromType() or
/// ArrowSchemaDeepCopy().
ArrowErrorCode ArrowSchemaAllocateChildren(struct ArrowSchema* schema,
                                           int64_t n_children);

/// \brief Allocate the schema->dictionary member
///
/// schema must have been allocated using ArrowSchemaInitFromType() or
/// ArrowSchemaDeepCopy().
ArrowErrorCode ArrowSchemaAllocateDictionary(struct ArrowSchema* schema);

/// @}

/// \defgroup nanoarrow-metadata Create, read, and modify schema metadata
///
/// @{

/// \brief Reader for key/value pairs in schema metadata
///
/// The ArrowMetadataReader does not own any data and is only valid
/// for the lifetime of the underlying metadata pointer.
struct ArrowMetadataReader {
  /// \brief A metadata string from a schema->metadata field.
  const char* metadata;

  /// \brief The current offset into the metadata string
  int64_t offset;

  /// \brief The number of remaining keys
  int32_t remaining_keys;
};

/// \brief Initialize an ArrowMetadataReader
ArrowErrorCode ArrowMetadataReaderInit(struct ArrowMetadataReader* reader,
                                       const char* metadata);

/// \brief Read the next key/value pair from an ArrowMetadataReader
ArrowErrorCode ArrowMetadataReaderRead(struct ArrowMetadataReader* reader,
                                       struct ArrowStringView* key_out,
                                       struct ArrowStringView* value_out);

/// \brief The number of bytes in in a key/value metadata string
int64_t ArrowMetadataSizeOf(const char* metadata);

/// \brief Check for a key in schema metadata
char ArrowMetadataHasKey(const char* metadata, struct ArrowStringView key);

/// \brief Extract a value from schema metadata
///
/// If key does not exist in metadata, value_out is unmodified
ArrowErrorCode ArrowMetadataGetValue(const char* metadata, struct ArrowStringView key,
                                     struct ArrowStringView* value_out);

/// \brief Initialize a builder for schema metadata from key/value pairs
///
/// metadata can be an existing metadata string or NULL to initialize
/// an empty metadata string.
ArrowErrorCode ArrowMetadataBuilderInit(struct ArrowBuffer* buffer, const char* metadata);

/// \brief Append a key/value pair to a buffer containing serialized metadata
ArrowErrorCode ArrowMetadataBuilderAppend(struct ArrowBuffer* buffer,
                                          struct ArrowStringView key,
                                          struct ArrowStringView value);

/// \brief Set a key/value pair to a buffer containing serialized metadata
///
/// Ensures that the only entry for key in the metadata is set to value.
/// This function maintains the existing position of (the first instance of)
/// key if present in the data.
ArrowErrorCode ArrowMetadataBuilderSet(struct ArrowBuffer* buffer,
                                       struct ArrowStringView key,
                                       struct ArrowStringView value);

/// \brief Remove a key from a buffer containing serialized metadata
ArrowErrorCode ArrowMetadataBuilderRemove(struct ArrowBuffer* buffer,
                                          struct ArrowStringView key);

/// @}

/// \defgroup nanoarrow-schema-view Reading schemas
///
/// @{

/// \brief A non-owning view of a parsed ArrowSchema
///
/// Contains more readily extractable values than a raw ArrowSchema.
/// Clients can stack or statically allocate this structure but are
/// encouraged to use the provided getters to ensure forward
/// compatibility.
struct ArrowSchemaView {
  /// \brief A pointer to the schema represented by this view
  struct ArrowSchema* schema;

  /// \brief The data type represented by the schema
  ///
  /// This value may be NANOARROW_TYPE_DICTIONARY if the schema has a
  /// non-null dictionary member; datetime types are valid values.
  /// This value will never be NANOARROW_TYPE_EXTENSION (see
  /// extension_name and/or extension_metadata to check for
  /// an extension type).
  enum ArrowType type;

  /// \brief The storage data type represented by the schema
  ///
  /// This value will never be NANOARROW_TYPE_DICTIONARY, NANOARROW_TYPE_EXTENSION
  /// or any datetime type. This value represents only the type required to
  /// interpret the buffers in the array.
  enum ArrowType storage_type;

  /// \brief The storage layout represented by the schema
  struct ArrowLayout layout;

  /// \brief The extension type name if it exists
  ///
  /// If the ARROW:extension:name key is present in schema.metadata,
  /// extension_name.data will be non-NULL.
  struct ArrowStringView extension_name;

  /// \brief The extension type metadata if it exists
  ///
  /// If the ARROW:extension:metadata key is present in schema.metadata,
  /// extension_metadata.data will be non-NULL.
  struct ArrowStringView extension_metadata;

  /// \brief Format fixed size parameter
  ///
  /// This value is set when parsing a fixed-size binary or fixed-size
  /// list schema; this value is undefined for other types. For a
  /// fixed-size binary schema this value is in bytes; for a fixed-size
  /// list schema this value refers to the number of child elements for
  /// each element of the parent.
  int32_t fixed_size;

  /// \brief Decimal bitwidth
  ///
  /// This value is set when parsing a decimal type schema;
  /// this value is undefined for other types.
  int32_t decimal_bitwidth;

  /// \brief Decimal precision
  ///
  /// This value is set when parsing a decimal type schema;
  /// this value is undefined for other types.
  int32_t decimal_precision;

  /// \brief Decimal scale
  ///
  /// This value is set when parsing a decimal type schema;
  /// this value is undefined for other types.
  int32_t decimal_scale;

  /// \brief Format time unit parameter
  ///
  /// This value is set when parsing a date/time type. The value is
  /// undefined for other types.
  enum ArrowTimeUnit time_unit;

  /// \brief Format timezone parameter
  ///
  /// This value is set when parsing a timestamp type and represents
  /// the timezone format parameter. This value points to
  /// data within the schema and is undefined for other types.
  const char* timezone;

  /// \brief Union type ids parameter
  ///
  /// This value is set when parsing a union type and represents
  /// type ids parameter. This value points to
  /// data within the schema and is undefined for other types.
  const char* union_type_ids;
};

/// \brief Initialize an ArrowSchemaView
ArrowErrorCode ArrowSchemaViewInit(struct ArrowSchemaView* schema_view,
                                   struct ArrowSchema* schema, struct ArrowError* error);

/// @}

/// \defgroup nanoarrow-buffer Owning, growable buffers
///
/// @{

/// \brief Initialize an ArrowBuffer
///
/// Initialize a buffer with a NULL, zero-size buffer using the default
/// buffer allocator.
static inline void ArrowBufferInit(struct ArrowBuffer* buffer);

/// \brief Set a newly-initialized buffer's allocator
///
/// Returns EINVAL if the buffer has already been allocated.
static inline ArrowErrorCode ArrowBufferSetAllocator(
    struct ArrowBuffer* buffer, struct ArrowBufferAllocator allocator);

/// \brief Reset an ArrowBuffer
///
/// Releases the buffer using the allocator's free method if
/// the buffer's data member is non-null, sets the data member
/// to NULL, and sets the buffer's size and capacity to 0.
static inline void ArrowBufferReset(struct ArrowBuffer* buffer);

/// \brief Move an ArrowBuffer
///
/// Transfers the buffer data and lifecycle management to another
/// address and resets buffer.
static inline void ArrowBufferMove(struct ArrowBuffer* src, struct ArrowBuffer* dst);

/// \brief Grow or shrink a buffer to a given capacity
///
/// When shrinking the capacity of the buffer, the buffer is only reallocated
/// if shrink_to_fit is non-zero. Calling ArrowBufferResize() does not
/// adjust the buffer's size member except to ensure that the invariant
/// capacity >= size remains true.
static inline ArrowErrorCode ArrowBufferResize(struct ArrowBuffer* buffer,
                                               int64_t new_capacity_bytes,
                                               char shrink_to_fit);

/// \brief Ensure a buffer has at least a given additional capacity
///
/// Ensures that the buffer has space to append at least
/// additional_size_bytes, overallocating when required.
static inline ArrowErrorCode ArrowBufferReserve(struct ArrowBuffer* buffer,
                                                int64_t additional_size_bytes);

/// \brief Write data to buffer and increment the buffer size
///
/// This function does not check that buffer has the required capacity
static inline void ArrowBufferAppendUnsafe(struct ArrowBuffer* buffer, const void* data,
                                           int64_t size_bytes);

/// \brief Write data to buffer and increment the buffer size
///
/// This function writes and ensures that the buffer has the required capacity,
/// possibly by reallocating the buffer. Like ArrowBufferReserve, this will
/// overallocate when reallocation is required.
static inline ArrowErrorCode ArrowBufferAppend(struct ArrowBuffer* buffer,
                                               const void* data, int64_t size_bytes);

/// \brief Write fill to buffer and increment the buffer size
///
/// This function writes the specified number of fill bytes and
/// ensures that the buffer has the required capacity,
static inline ArrowErrorCode ArrowBufferAppendFill(struct ArrowBuffer* buffer,
                                                   uint8_t value, int64_t size_bytes);

/// \brief Write an 8-bit integer to a buffer
static inline ArrowErrorCode ArrowBufferAppendInt8(struct ArrowBuffer* buffer,
                                                   int8_t value);

/// \brief Write an unsigned 8-bit integer to a buffer
static inline ArrowErrorCode ArrowBufferAppendUInt8(struct ArrowBuffer* buffer,
                                                    uint8_t value);

/// \brief Write a 16-bit integer to a buffer
static inline ArrowErrorCode ArrowBufferAppendInt16(struct ArrowBuffer* buffer,
                                                    int16_t value);

/// \brief Write an unsigned 16-bit integer to a buffer
static inline ArrowErrorCode ArrowBufferAppendUInt16(struct ArrowBuffer* buffer,
                                                     uint16_t value);

/// \brief Write a 32-bit integer to a buffer
static inline ArrowErrorCode ArrowBufferAppendInt32(struct ArrowBuffer* buffer,
                                                    int32_t value);

/// \brief Write an unsigned 32-bit integer to a buffer
static inline ArrowErrorCode ArrowBufferAppendUInt32(struct ArrowBuffer* buffer,
                                                     uint32_t value);

/// \brief Write a 64-bit integer to a buffer
static inline ArrowErrorCode ArrowBufferAppendInt64(struct ArrowBuffer* buffer,
                                                    int64_t value);

/// \brief Write an unsigned 64-bit integer to a buffer
static inline ArrowErrorCode ArrowBufferAppendUInt64(struct ArrowBuffer* buffer,
                                                     uint64_t value);

/// \brief Write a double to a buffer
static inline ArrowErrorCode ArrowBufferAppendDouble(struct ArrowBuffer* buffer,
                                                     double value);

/// \brief Write a float to a buffer
static inline ArrowErrorCode ArrowBufferAppendFloat(struct ArrowBuffer* buffer,
                                                    float value);

/// \brief Write an ArrowStringView to a buffer
static inline ArrowErrorCode ArrowBufferAppendStringView(struct ArrowBuffer* buffer,
                                                         struct ArrowStringView value);

/// \brief Write an ArrowBufferView to a buffer
static inline ArrowErrorCode ArrowBufferAppendBufferView(struct ArrowBuffer* buffer,
                                                         struct ArrowBufferView value);

/// @}

/// \defgroup nanoarrow-bitmap Bitmap utilities
///
/// @{

/// \brief Extract a boolean value from a bitmap
static inline int8_t ArrowBitGet(const uint8_t* bits, int64_t i);

/// \brief Set a boolean value to a bitmap to true
static inline void ArrowBitSet(uint8_t* bits, int64_t i);

/// \brief Set a boolean value to a bitmap to false
static inline void ArrowBitClear(uint8_t* bits, int64_t i);

/// \brief Set a boolean value to a bitmap
static inline void ArrowBitSetTo(uint8_t* bits, int64_t i, uint8_t value);

/// \brief Set a boolean value to a range in a bitmap
static inline void ArrowBitsSetTo(uint8_t* bits, int64_t start_offset, int64_t length,
                                  uint8_t bits_are_set);

/// \brief Count true values in a bitmap
static inline int64_t ArrowBitCountSet(const uint8_t* bits, int64_t i_from, int64_t i_to);

/// \brief Initialize an ArrowBitmap
///
/// Initialize the builder's buffer, empty its cache, and reset the size to zero
static inline void ArrowBitmapInit(struct ArrowBitmap* bitmap);

/// \brief Move an ArrowBitmap
///
/// Transfers the underlying buffer data and lifecycle management to another
/// address and resets the bitmap.
static inline void ArrowBitmapMove(struct ArrowBitmap* src, struct ArrowBitmap* dst);

/// \brief Ensure a bitmap builder has at least a given additional capacity
///
/// Ensures that the buffer has space to append at least
/// additional_size_bits, overallocating when required.
static inline ArrowErrorCode ArrowBitmapReserve(struct ArrowBitmap* bitmap,
                                                int64_t additional_size_bits);

/// \brief Grow or shrink a bitmap to a given capacity
///
/// When shrinking the capacity of the bitmap, the bitmap is only reallocated
/// if shrink_to_fit is non-zero. Calling ArrowBitmapResize() does not
/// adjust the buffer's size member except when shrinking new_capacity_bits
/// to a value less than the current number of bits in the bitmap.
static inline ArrowErrorCode ArrowBitmapResize(struct ArrowBitmap* bitmap,
                                               int64_t new_capacity_bits,
                                               char shrink_to_fit);

/// \brief Reserve space for and append zero or more of the same boolean value to a bitmap
static inline ArrowErrorCode ArrowBitmapAppend(struct ArrowBitmap* bitmap,
                                               uint8_t bits_are_set, int64_t length);

/// \brief Append zero or more of the same boolean value to a bitmap
static inline void ArrowBitmapAppendUnsafe(struct ArrowBitmap* bitmap,
                                           uint8_t bits_are_set, int64_t length);

/// \brief Append boolean values encoded as int8_t to a bitmap
///
/// The values must all be 0 or 1.
static inline void ArrowBitmapAppendInt8Unsafe(struct ArrowBitmap* bitmap,
                                               const int8_t* values, int64_t n_values);

/// \brief Append boolean values encoded as int32_t to a bitmap
///
/// The values must all be 0 or 1.
static inline void ArrowBitmapAppendInt32Unsafe(struct ArrowBitmap* bitmap,
                                                const int32_t* values, int64_t n_values);

/// \brief Reset a bitmap builder
///
/// Releases any memory held by buffer, empties the cache, and resets the size to zero
static inline void ArrowBitmapReset(struct ArrowBitmap* bitmap);

/// @}

/// \defgroup nanoarrow-array Creating arrays
///
/// These functions allocate, copy, and destroy ArrowArray structures.
/// Once an ArrowArray has been initialized via ArrowArrayInitFromType()
/// or ArrowArrayInitFromSchema(), the caller is responsible for releasing
/// it using the embedded release callback.
///
/// @{

/// \brief Initialize the fields of an array
///
/// Initializes the fields and release callback of array. Caller
/// is responsible for calling the array->release callback if
/// NANOARROW_OK is returned.
ArrowErrorCode ArrowArrayInitFromType(struct ArrowArray* array,
                                      enum ArrowType storage_type);

/// \brief Initialize the contents of an ArrowArray from an ArrowSchema
///
/// Caller is responsible for calling the array->release callback if
/// NANOARROW_OK is returned.
ArrowErrorCode ArrowArrayInitFromSchema(struct ArrowArray* array,
                                        struct ArrowSchema* schema,
                                        struct ArrowError* error);

/// \brief Initialize the contents of an ArrowArray from an ArrowArrayView
///
/// Caller is responsible for calling the array->release callback if
/// NANOARROW_OK is returned.
ArrowErrorCode ArrowArrayInitFromArrayView(struct ArrowArray* array,
                                           struct ArrowArrayView* array_view,
                                           struct ArrowError* error);

/// \brief Allocate the array->children array
///
/// Includes the memory for each child struct ArrowArray,
/// whose members are marked as released and may be subsequently initialized
/// with ArrowArrayInitFromType() or moved from an existing ArrowArray.
/// schema must have been allocated using ArrowArrayInitFromType().
ArrowErrorCode ArrowArrayAllocateChildren(struct ArrowArray* array, int64_t n_children);

/// \brief Allocate the array->dictionary member
///
/// Includes the memory for the struct ArrowArray, whose contents
/// is marked as released and may be subsequently initialized
/// with ArrowArrayInitFromType() or moved from an existing ArrowArray.
/// array must have been allocated using ArrowArrayInitFromType()
ArrowErrorCode ArrowArrayAllocateDictionary(struct ArrowArray* array);

/// \brief Set the validity bitmap of an ArrowArray
///
/// array must have been allocated using ArrowArrayInitFromType()
void ArrowArraySetValidityBitmap(struct ArrowArray* array, struct ArrowBitmap* bitmap);

/// \brief Set a buffer of an ArrowArray
///
/// array must have been allocated using ArrowArrayInitFromType()
ArrowErrorCode ArrowArraySetBuffer(struct ArrowArray* array, int64_t i,
                                   struct ArrowBuffer* buffer);

/// \brief Get the validity bitmap of an ArrowArray
///
/// array must have been allocated using ArrowArrayInitFromType()
static inline struct ArrowBitmap* ArrowArrayValidityBitmap(struct ArrowArray* array);

/// \brief Get a buffer of an ArrowArray
///
/// array must have been allocated using ArrowArrayInitFromType()
static inline struct ArrowBuffer* ArrowArrayBuffer(struct ArrowArray* array, int64_t i);

/// \brief Start element-wise appending to an ArrowArray
///
/// Initializes any values needed to use ArrowArrayAppend*() functions.
/// All element-wise appenders append by value and return EINVAL if the exact value
/// cannot be represented by the underlying storage type.
/// array must have been allocated using ArrowArrayInitFromType()
static inline ArrowErrorCode ArrowArrayStartAppending(struct ArrowArray* array);

/// \brief Reserve space for future appends
///
/// For buffer sizes that can be calculated (i.e., not string data buffers or
/// child array sizes for non-fixed-size arrays), recursively reserve space for
/// additional elements. This is useful for reducing the number of reallocations
/// that occur using the item-wise appenders.
ArrowErrorCode ArrowArrayReserve(struct ArrowArray* array,
                                 int64_t additional_size_elements);

/// \brief Append a null value to an array
static inline ArrowErrorCode ArrowArrayAppendNull(struct ArrowArray* array, int64_t n);

/// \brief Append an empty, non-null value to an array
static inline ArrowErrorCode ArrowArrayAppendEmpty(struct ArrowArray* array, int64_t n);

/// \brief Append a signed integer value to an array
///
/// Returns NANOARROW_OK if value can be exactly represented by
/// the underlying storage type or EINVAL otherwise (e.g., value
/// is outside the valid array range).
static inline ArrowErrorCode ArrowArrayAppendInt(struct ArrowArray* array, int64_t value);

/// \brief Append an unsigned integer value to an array
///
/// Returns NANOARROW_OK if value can be exactly represented by
/// the underlying storage type or EINVAL otherwise (e.g., value
/// is outside the valid array range).
static inline ArrowErrorCode ArrowArrayAppendUInt(struct ArrowArray* array,
                                                  uint64_t value);

/// \brief Append a double value to an array
///
/// Returns NANOARROW_OK if value can be exactly represented by
/// the underlying storage type or EINVAL otherwise (e.g., value
/// is outside the valid array range or there is an attempt to append
/// a non-integer to an array with an integer storage type).
static inline ArrowErrorCode ArrowArrayAppendDouble(struct ArrowArray* array,
                                                    double value);

/// \brief Append a string of bytes to an array
///
/// Returns NANOARROW_OK if value can be exactly represented by
/// the underlying storage type or EINVAL otherwise (e.g.,
/// the underlying array is not a binary, string, large binary, large string,
/// or fixed-size binary array, or value is the wrong size for a fixed-size
/// binary array).
static inline ArrowErrorCode ArrowArrayAppendBytes(struct ArrowArray* array,
                                                   struct ArrowBufferView value);

/// \brief Append a string value to an array
///
/// Returns NANOARROW_OK if value can be exactly represented by
/// the underlying storage type or EINVAL otherwise (e.g.,
/// the underlying array is not a string or large string array).
static inline ArrowErrorCode ArrowArrayAppendString(struct ArrowArray* array,
                                                    struct ArrowStringView value);

/// \brief Append a Interval to an array
///
/// Returns NANOARROW_OK if value can be exactly represented by
/// the underlying storage type or EINVAL otherwise.
static inline ArrowErrorCode ArrowArrayAppendInterval(struct ArrowArray* array,
                                                      struct ArrowInterval* value);

/// \brief Append a decimal value to an array
///
/// Returns NANOARROW_OK if array is a decimal array with the appropriate
/// bitwidth or EINVAL otherwise.
static inline ArrowErrorCode ArrowArrayAppendDecimal(struct ArrowArray* array,
                                                     struct ArrowDecimal* value);

/// \brief Finish a nested array element
///
/// Appends a non-null element to the array based on the first child's current
/// length. Returns NANOARROW_OK if the item was successfully added or EINVAL
/// if the underlying storage type is not a struct, list, large list, or fixed-size
/// list, or if there was an attempt to add a struct or fixed-size list element where the
/// length of the child array(s) did not match the expected length.
static inline ArrowErrorCode ArrowArrayFinishElement(struct ArrowArray* array);

/// \brief Finish a union array element
///
/// Appends an element to the union type ids buffer and increments array->length.
/// For sparse unions, up to one element is added to non type-id children. Returns
/// EINVAL if the underlying storage type is not a union, if type_id is not valid,
/// or if child sizes after appending are inconsistent.
static inline ArrowErrorCode ArrowArrayFinishUnionElement(struct ArrowArray* array,
                                                          int8_t type_id);

/// \brief Shrink buffer capacity to the size required
///
/// Also applies shrinking to any child arrays. array must have been allocated using
/// ArrowArrayInitFromType
static inline ArrowErrorCode ArrowArrayShrinkToFit(struct ArrowArray* array);

/// \brief Finish building an ArrowArray
///
/// Flushes any pointers from internal buffers that may have been reallocated
/// into array->buffers and checks the actual size of the buffers
/// against the expected size based on the final length.
/// array must have been allocated using ArrowArrayInitFromType()
ArrowErrorCode ArrowArrayFinishBuildingDefault(struct ArrowArray* array,
                                               struct ArrowError* error);

/// \brief Finish building an ArrowArray with explicit validation
///
/// Finish building with an explicit validation level. This could perform less validation
/// (i.e. NANOARROW_VALIDATION_LEVEL_NONE or NANOARROW_VALIDATION_LEVEL_MINIMAL) if CPU
/// buffer data access is not possible or more validation (i.e.,
/// NANOARROW_VALIDATION_LEVEL_FULL) if buffer content was obtained from an untrusted or
/// corruptible source.
ArrowErrorCode ArrowArrayFinishBuilding(struct ArrowArray* array,
                                        enum ArrowValidationLevel validation_level,
                                        struct ArrowError* error);

/// @}

/// \defgroup nanoarrow-array-view Reading arrays
///
/// These functions read and validate the contents ArrowArray structures.
///
/// @{

/// \brief Initialize the contents of an ArrowArrayView
void ArrowArrayViewInitFromType(struct ArrowArrayView* array_view,
                                enum ArrowType storage_type);

/// \brief Move an ArrowArrayView
///
/// Transfers the ArrowArrayView data and lifecycle management to another
/// address and resets the contents of src.
static inline void ArrowArrayViewMove(struct ArrowArrayView* src,
                                      struct ArrowArrayView* dst);

/// \brief Initialize the contents of an ArrowArrayView from an ArrowSchema
ArrowErrorCode ArrowArrayViewInitFromSchema(struct ArrowArrayView* array_view,
                                            struct ArrowSchema* schema,
                                            struct ArrowError* error);

/// \brief Allocate the array_view->children array
///
/// Includes the memory for each child struct ArrowArrayView
ArrowErrorCode ArrowArrayViewAllocateChildren(struct ArrowArrayView* array_view,
                                              int64_t n_children);

/// \brief Allocate array_view->dictionary
ArrowErrorCode ArrowArrayViewAllocateDictionary(struct ArrowArrayView* array_view);

/// \brief Set data-independent buffer sizes from length
void ArrowArrayViewSetLength(struct ArrowArrayView* array_view, int64_t length);

/// \brief Set buffer sizes and data pointers from an ArrowArray
ArrowErrorCode ArrowArrayViewSetArray(struct ArrowArrayView* array_view,
                                      struct ArrowArray* array, struct ArrowError* error);

/// \brief Set buffer sizes and data pointers from an ArrowArray except for those
/// that require dereferencing buffer content.
ArrowErrorCode ArrowArrayViewSetArrayMinimal(struct ArrowArrayView* array_view,
                                             struct ArrowArray* array,
                                             struct ArrowError* error);

/// \brief Performs checks on the content of an ArrowArrayView
///
/// If using ArrowArrayViewSetArray() to back array_view with an ArrowArray,
/// the buffer sizes and some content (fist and last offset) have already
/// been validated at the "default" level. If setting the buffer pointers
/// and sizes otherwise, you may wish to perform checks at a different level. See
/// documentation for ArrowValidationLevel for the details of checks performed
/// at each level.
ArrowErrorCode ArrowArrayViewValidate(struct ArrowArrayView* array_view,
                                      enum ArrowValidationLevel validation_level,
                                      struct ArrowError* error);

/// \brief Reset the contents of an ArrowArrayView and frees resources
void ArrowArrayViewReset(struct ArrowArrayView* array_view);

/// \brief Check for a null element in an ArrowArrayView
static inline int8_t ArrowArrayViewIsNull(struct ArrowArrayView* array_view, int64_t i);

/// \brief Get the type id of a union array element
static inline int8_t ArrowArrayViewUnionTypeId(struct ArrowArrayView* array_view,
                                               int64_t i);

/// \brief Get the child index of a union array element
static inline int8_t ArrowArrayViewUnionChildIndex(struct ArrowArrayView* array_view,
                                                   int64_t i);

/// \brief Get the index to use into the relevant union child array
static inline int64_t ArrowArrayViewUnionChildOffset(struct ArrowArrayView* array_view,
                                                     int64_t i);

/// \brief Get an element in an ArrowArrayView as an integer
///
/// This function does not check for null values, that values are actually integers, or
/// that values are within a valid range for an int64.
static inline int64_t ArrowArrayViewGetIntUnsafe(struct ArrowArrayView* array_view,
                                                 int64_t i);

/// \brief Get an element in an ArrowArrayView as an unsigned integer
///
/// This function does not check for null values, that values are actually integers, or
/// that values are within a valid range for a uint64.
static inline uint64_t ArrowArrayViewGetUIntUnsafe(struct ArrowArrayView* array_view,
                                                   int64_t i);

/// \brief Get an element in an ArrowArrayView as a double
///
/// This function does not check for null values, or
/// that values are within a valid range for a double.
static inline double ArrowArrayViewGetDoubleUnsafe(struct ArrowArrayView* array_view,
                                                   int64_t i);

/// \brief Get an element in an ArrowArrayView as an ArrowStringView
///
/// This function does not check for null values.
static inline struct ArrowStringView ArrowArrayViewGetStringUnsafe(
    struct ArrowArrayView* array_view, int64_t i);

/// \brief Get an element in an ArrowArrayView as an ArrowBufferView
///
/// This function does not check for null values.
static inline struct ArrowBufferView ArrowArrayViewGetBytesUnsafe(
    struct ArrowArrayView* array_view, int64_t i);

/// \brief Get an element in an ArrowArrayView as an ArrowDecimal
///
/// This function does not check for null values. The out parameter must
/// be initialized with ArrowDecimalInit() with the proper parameters for this
/// type before calling this for the first time.
static inline void ArrowArrayViewGetDecimalUnsafe(struct ArrowArrayView* array_view,
                                                  int64_t i, struct ArrowDecimal* out);

/// @}

/// \defgroup nanoarrow-basic-array-stream Basic ArrowArrayStream implementation
///
/// An implementation of an ArrowArrayStream based on a collection of
/// zero or more previously-existing ArrowArray objects. Users should
/// initialize and/or validate the contents before transferring the
/// responsibility of the ArrowArrayStream elsewhere.
///
/// @{

/// \brief Initialize an ArrowArrayStream backed by this implementation
///
/// This function moves the ownership of schema to the array_stream. If
/// this function returns NANOARROW_OK, the caller is responsible for
/// releasing the ArrowArrayStream.
ArrowErrorCode ArrowBasicArrayStreamInit(struct ArrowArrayStream* array_stream,
                                         struct ArrowSchema* schema, int64_t n_arrays);

/// \brief Set the ith ArrowArray in this ArrowArrayStream.
///
/// array_stream must have been initialized with ArrowBasicArrayStreamInit().
/// This function move the ownership of array to the array_stream. i must
/// be greater than zero and less than the value of n_arrays passed in
/// ArrowBasicArrayStreamInit(). Callers are not required to fill all
/// n_arrays members (i.e., n_arrays is a maximum bound).
void ArrowBasicArrayStreamSetArray(struct ArrowArrayStream* array_stream, int64_t i,
                                   struct ArrowArray* array);

/// \brief Validate the contents of this ArrowArrayStream
///
/// array_stream must have been initialized with ArrowBasicArrayStreamInit().
/// This function uses ArrowArrayStreamInitFromSchema() and ArrowArrayStreamSetArray()
/// to validate the contents of the arrays.
ArrowErrorCode ArrowBasicArrayStreamValidate(struct ArrowArrayStream* array_stream,
                                             struct ArrowError* error);

/// @}

// Inline function definitions



#ifdef __cplusplus
}
#endif

#endif
// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#ifndef NANOARROW_BUFFER_INLINE_H_INCLUDED
#define NANOARROW_BUFFER_INLINE_H_INCLUDED

#include <errno.h>
#include <stdint.h>
#include <string.h>



#ifdef __cplusplus
extern "C" {
#endif

static inline int64_t _ArrowGrowByFactor(int64_t current_capacity, int64_t new_capacity) {
  int64_t doubled_capacity = current_capacity * 2;
  if (doubled_capacity > new_capacity) {
    return doubled_capacity;
  } else {
    return new_capacity;
  }
}

static inline void ArrowBufferInit(struct ArrowBuffer* buffer) {
  buffer->data = NULL;
  buffer->size_bytes = 0;
  buffer->capacity_bytes = 0;
  buffer->allocator = ArrowBufferAllocatorDefault();
}

static inline ArrowErrorCode ArrowBufferSetAllocator(
    struct ArrowBuffer* buffer, struct ArrowBufferAllocator allocator) {
  if (buffer->data == NULL) {
    buffer->allocator = allocator;
    return NANOARROW_OK;
  } else {
    return EINVAL;
  }
}

static inline void ArrowBufferReset(struct ArrowBuffer* buffer) {
  if (buffer->data != NULL) {
    buffer->allocator.free(&buffer->allocator, (uint8_t*)buffer->data,
                           buffer->capacity_bytes);
    buffer->data = NULL;
  }

  buffer->capacity_bytes = 0;
  buffer->size_bytes = 0;
}

static inline void ArrowBufferMove(struct ArrowBuffer* src, struct ArrowBuffer* dst) {
  memcpy(dst, src, sizeof(struct ArrowBuffer));
  src->data = NULL;
  ArrowBufferReset(src);
}

static inline ArrowErrorCode ArrowBufferResize(struct ArrowBuffer* buffer,
                                               int64_t new_capacity_bytes,
                                               char shrink_to_fit) {
  if (new_capacity_bytes < 0) {
    return EINVAL;
  }

  if (new_capacity_bytes > buffer->capacity_bytes || shrink_to_fit) {
    buffer->data = buffer->allocator.reallocate(
        &buffer->allocator, buffer->data, buffer->capacity_bytes, new_capacity_bytes);
    if (buffer->data == NULL && new_capacity_bytes > 0) {
      buffer->capacity_bytes = 0;
      buffer->size_bytes = 0;
      return ENOMEM;
    }

    buffer->capacity_bytes = new_capacity_bytes;
  }

  // Ensures that when shrinking that size <= capacity
  if (new_capacity_bytes < buffer->size_bytes) {
    buffer->size_bytes = new_capacity_bytes;
  }

  return NANOARROW_OK;
}

static inline ArrowErrorCode ArrowBufferReserve(struct ArrowBuffer* buffer,
                                                int64_t additional_size_bytes) {
  int64_t min_capacity_bytes = buffer->size_bytes + additional_size_bytes;
  if (min_capacity_bytes <= buffer->capacity_bytes) {
    return NANOARROW_OK;
  }

  return ArrowBufferResize(
      buffer, _ArrowGrowByFactor(buffer->capacity_bytes, min_capacity_bytes), 0);
}

static inline void ArrowBufferAppendUnsafe(struct ArrowBuffer* buffer, const void* data,
                                           int64_t size_bytes) {
  if (size_bytes > 0) {
    memcpy(buffer->data + buffer->size_bytes, data, size_bytes);
    buffer->size_bytes += size_bytes;
  }
}

static inline ArrowErrorCode ArrowBufferAppend(struct ArrowBuffer* buffer,
                                               const void* data, int64_t size_bytes) {
  NANOARROW_RETURN_NOT_OK(ArrowBufferReserve(buffer, size_bytes));

  ArrowBufferAppendUnsafe(buffer, data, size_bytes);
  return NANOARROW_OK;
}

static inline ArrowErrorCode ArrowBufferAppendInt8(struct ArrowBuffer* buffer,
                                                   int8_t value) {
  return ArrowBufferAppend(buffer, &value, sizeof(int8_t));
}

static inline ArrowErrorCode ArrowBufferAppendUInt8(struct ArrowBuffer* buffer,
                                                    uint8_t value) {
  return ArrowBufferAppend(buffer, &value, sizeof(uint8_t));
}

static inline ArrowErrorCode ArrowBufferAppendInt16(struct ArrowBuffer* buffer,
                                                    int16_t value) {
  return ArrowBufferAppend(buffer, &value, sizeof(int16_t));
}

static inline ArrowErrorCode ArrowBufferAppendUInt16(struct ArrowBuffer* buffer,
                                                     uint16_t value) {
  return ArrowBufferAppend(buffer, &value, sizeof(uint16_t));
}

static inline ArrowErrorCode ArrowBufferAppendInt32(struct ArrowBuffer* buffer,
                                                    int32_t value) {
  return ArrowBufferAppend(buffer, &value, sizeof(int32_t));
}

static inline ArrowErrorCode ArrowBufferAppendUInt32(struct ArrowBuffer* buffer,
                                                     uint32_t value) {
  return ArrowBufferAppend(buffer, &value, sizeof(uint32_t));
}

static inline ArrowErrorCode ArrowBufferAppendInt64(struct ArrowBuffer* buffer,
                                                    int64_t value) {
  return ArrowBufferAppend(buffer, &value, sizeof(int64_t));
}

static inline ArrowErrorCode ArrowBufferAppendUInt64(struct ArrowBuffer* buffer,
                                                     uint64_t value) {
  return ArrowBufferAppend(buffer, &value, sizeof(uint64_t));
}

static inline ArrowErrorCode ArrowBufferAppendDouble(struct ArrowBuffer* buffer,
                                                     double value) {
  return ArrowBufferAppend(buffer, &value, sizeof(double));
}

static inline ArrowErrorCode ArrowBufferAppendFloat(struct ArrowBuffer* buffer,
                                                    float value) {
  return ArrowBufferAppend(buffer, &value, sizeof(float));
}

static inline ArrowErrorCode ArrowBufferAppendStringView(struct ArrowBuffer* buffer,
                                                         struct ArrowStringView value) {
  return ArrowBufferAppend(buffer, value.data, value.size_bytes);
}

static inline ArrowErrorCode ArrowBufferAppendBufferView(struct ArrowBuffer* buffer,
                                                         struct ArrowBufferView value) {
  return ArrowBufferAppend(buffer, value.data.data, value.size_bytes);
}

static inline ArrowErrorCode ArrowBufferAppendFill(struct ArrowBuffer* buffer,
                                                   uint8_t value, int64_t size_bytes) {
  NANOARROW_RETURN_NOT_OK(ArrowBufferReserve(buffer, size_bytes));

  memset(buffer->data + buffer->size_bytes, value, size_bytes);
  buffer->size_bytes += size_bytes;
  return NANOARROW_OK;
}

static const uint8_t _ArrowkBitmask[] = {1, 2, 4, 8, 16, 32, 64, 128};
static const uint8_t _ArrowkFlippedBitmask[] = {254, 253, 251, 247, 239, 223, 191, 127};
static const uint8_t _ArrowkPrecedingBitmask[] = {0, 1, 3, 7, 15, 31, 63, 127};
static const uint8_t _ArrowkTrailingBitmask[] = {255, 254, 252, 248, 240, 224, 192, 128};

static const uint8_t _ArrowkBytePopcount[] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3,
    4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4,
    4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4,
    5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5,
    4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2,
    3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5,
    5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4,
    5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};

static inline int64_t _ArrowRoundUpToMultipleOf8(int64_t value) {
  return (value + 7) & ~((int64_t)7);
}

static inline int64_t _ArrowRoundDownToMultipleOf8(int64_t value) {
  return (value / 8) * 8;
}

static inline int64_t _ArrowBytesForBits(int64_t bits) {
  return (bits >> 3) + ((bits & 7) != 0);
}

static inline void _ArrowBitmapPackInt8(const int8_t* values, uint8_t* out) {
  *out = (values[0] | values[1] << 1 | values[2] << 2 | values[3] << 3 | values[4] << 4 |
          values[5] << 5 | values[6] << 6 | values[7] << 7);
}

static inline void _ArrowBitmapPackInt32(const int32_t* values, uint8_t* out) {
  *out = (values[0] | values[1] << 1 | values[2] << 2 | values[3] << 3 | values[4] << 4 |
          values[5] << 5 | values[6] << 6 | values[7] << 7);
}

static inline int8_t ArrowBitGet(const uint8_t* bits, int64_t i) {
  return (bits[i >> 3] >> (i & 0x07)) & 1;
}

static inline void ArrowBitSet(uint8_t* bits, int64_t i) {
  bits[i / 8] |= _ArrowkBitmask[i % 8];
}

static inline void ArrowBitClear(uint8_t* bits, int64_t i) {
  bits[i / 8] &= _ArrowkFlippedBitmask[i % 8];
}

static inline void ArrowBitSetTo(uint8_t* bits, int64_t i, uint8_t bit_is_set) {
  bits[i / 8] ^=
      ((uint8_t)(-((uint8_t)(bit_is_set != 0)) ^ bits[i / 8])) & _ArrowkBitmask[i % 8];
}

static inline void ArrowBitsSetTo(uint8_t* bits, int64_t start_offset, int64_t length,
                                  uint8_t bits_are_set) {
  const int64_t i_begin = start_offset;
  const int64_t i_end = start_offset + length;
  const uint8_t fill_byte = (uint8_t)(-bits_are_set);

  const int64_t bytes_begin = i_begin / 8;
  const int64_t bytes_end = i_end / 8 + 1;

  const uint8_t first_byte_mask = _ArrowkPrecedingBitmask[i_begin % 8];
  const uint8_t last_byte_mask = _ArrowkTrailingBitmask[i_end % 8];

  if (bytes_end == bytes_begin + 1) {
    // set bits within a single byte
    const uint8_t only_byte_mask =
        i_end % 8 == 0 ? first_byte_mask : (uint8_t)(first_byte_mask | last_byte_mask);
    bits[bytes_begin] &= only_byte_mask;
    bits[bytes_begin] |= (uint8_t)(fill_byte & ~only_byte_mask);
    return;
  }

  // set/clear trailing bits of first byte
  bits[bytes_begin] &= first_byte_mask;
  bits[bytes_begin] |= (uint8_t)(fill_byte & ~first_byte_mask);

  if (bytes_end - bytes_begin > 2) {
    // set/clear whole bytes
    memset(bits + bytes_begin + 1, fill_byte, (size_t)(bytes_end - bytes_begin - 2));
  }

  if (i_end % 8 == 0) {
    return;
  }

  // set/clear leading bits of last byte
  bits[bytes_end - 1] &= last_byte_mask;
  bits[bytes_end - 1] |= (uint8_t)(fill_byte & ~last_byte_mask);
}

static inline int64_t ArrowBitCountSet(const uint8_t* bits, int64_t start_offset,
                                       int64_t length) {
  if (length == 0) {
    return 0;
  }

  const int64_t i_begin = start_offset;
  const int64_t i_end = start_offset + length;
  const int64_t i_last_valid = i_end - 1;

  const int64_t bytes_begin = i_begin / 8;
  const int64_t bytes_last_valid = i_last_valid / 8;

  if (bytes_begin == bytes_last_valid) {
    // count bits within a single byte
    const uint8_t first_byte_mask = _ArrowkPrecedingBitmask[i_end % 8];
    const uint8_t last_byte_mask = _ArrowkTrailingBitmask[i_begin % 8];

    const uint8_t only_byte_mask =
        i_end % 8 == 0 ? last_byte_mask : (uint8_t)(first_byte_mask & last_byte_mask);

    const uint8_t byte_masked = bits[bytes_begin] & only_byte_mask;
    return _ArrowkBytePopcount[byte_masked];
  }

  const uint8_t first_byte_mask = _ArrowkPrecedingBitmask[i_begin % 8];
  const uint8_t last_byte_mask = i_end % 8 == 0 ? 0 : _ArrowkTrailingBitmask[i_end % 8];
  int64_t count = 0;

  // first byte
  count += _ArrowkBytePopcount[bits[bytes_begin] & ~first_byte_mask];

  // middle bytes
  for (int64_t i = bytes_begin + 1; i < bytes_last_valid; i++) {
    count += _ArrowkBytePopcount[bits[i]];
  }

  // last byte
  count += _ArrowkBytePopcount[bits[bytes_last_valid] & ~last_byte_mask];

  return count;
}

static inline void ArrowBitmapInit(struct ArrowBitmap* bitmap) {
  ArrowBufferInit(&bitmap->buffer);
  bitmap->size_bits = 0;
}

static inline void ArrowBitmapMove(struct ArrowBitmap* src, struct ArrowBitmap* dst) {
  ArrowBufferMove(&src->buffer, &dst->buffer);
  dst->size_bits = src->size_bits;
  src->size_bits = 0;
}

static inline ArrowErrorCode ArrowBitmapReserve(struct ArrowBitmap* bitmap,
                                                int64_t additional_size_bits) {
  int64_t min_capacity_bits = bitmap->size_bits + additional_size_bits;
  if (min_capacity_bits <= (bitmap->buffer.capacity_bytes * 8)) {
    return NANOARROW_OK;
  }

  NANOARROW_RETURN_NOT_OK(
      ArrowBufferReserve(&bitmap->buffer, _ArrowBytesForBits(additional_size_bits)));

  bitmap->buffer.data[bitmap->buffer.capacity_bytes - 1] = 0;
  return NANOARROW_OK;
}

static inline ArrowErrorCode ArrowBitmapResize(struct ArrowBitmap* bitmap,
                                               int64_t new_capacity_bits,
                                               char shrink_to_fit) {
  if (new_capacity_bits < 0) {
    return EINVAL;
  }

  int64_t new_capacity_bytes = _ArrowBytesForBits(new_capacity_bits);
  NANOARROW_RETURN_NOT_OK(
      ArrowBufferResize(&bitmap->buffer, new_capacity_bytes, shrink_to_fit));

  if (new_capacity_bits < bitmap->size_bits) {
    bitmap->size_bits = new_capacity_bits;
  }

  return NANOARROW_OK;
}

static inline ArrowErrorCode ArrowBitmapAppend(struct ArrowBitmap* bitmap,
                                               uint8_t bits_are_set, int64_t length) {
  NANOARROW_RETURN_NOT_OK(ArrowBitmapReserve(bitmap, length));

  ArrowBitmapAppendUnsafe(bitmap, bits_are_set, length);
  return NANOARROW_OK;
}

static inline void ArrowBitmapAppendUnsafe(struct ArrowBitmap* bitmap,
                                           uint8_t bits_are_set, int64_t length) {
  ArrowBitsSetTo(bitmap->buffer.data, bitmap->size_bits, length, bits_are_set);
  bitmap->size_bits += length;
  bitmap->buffer.size_bytes = _ArrowBytesForBits(bitmap->size_bits);
}

static inline void ArrowBitmapAppendInt8Unsafe(struct ArrowBitmap* bitmap,
                                               const int8_t* values, int64_t n_values) {
  if (n_values == 0) {
    return;
  }

  const int8_t* values_cursor = values;
  int64_t n_remaining = n_values;
  int64_t out_i_cursor = bitmap->size_bits;
  uint8_t* out_cursor = bitmap->buffer.data + bitmap->size_bits / 8;

  // First byte
  if ((out_i_cursor % 8) != 0) {
    int64_t n_partial_bits = _ArrowRoundUpToMultipleOf8(out_i_cursor) - out_i_cursor;
    for (int i = 0; i < n_partial_bits; i++) {
      ArrowBitSetTo(bitmap->buffer.data, out_i_cursor++, values[i]);
    }

    out_cursor++;
    values_cursor += n_partial_bits;
    n_remaining -= n_partial_bits;
  }

  // Middle bytes
  int64_t n_full_bytes = n_remaining / 8;
  for (int64_t i = 0; i < n_full_bytes; i++) {
    _ArrowBitmapPackInt8(values_cursor, out_cursor);
    values_cursor += 8;
    out_cursor++;
  }

  // Last byte
  out_i_cursor += n_full_bytes * 8;
  n_remaining -= n_full_bytes * 8;
  if (n_remaining > 0) {
    // Zero out the last byte
    *out_cursor = 0x00;
    for (int i = 0; i < n_remaining; i++) {
      ArrowBitSetTo(bitmap->buffer.data, out_i_cursor++, values_cursor[i]);
    }
    out_cursor++;
  }

  bitmap->size_bits += n_values;
  bitmap->buffer.size_bytes = out_cursor - bitmap->buffer.data;
}

static inline void ArrowBitmapAppendInt32Unsafe(struct ArrowBitmap* bitmap,
                                                const int32_t* values, int64_t n_values) {
  if (n_values == 0) {
    return;
  }

  const int32_t* values_cursor = values;
  int64_t n_remaining = n_values;
  int64_t out_i_cursor = bitmap->size_bits;
  uint8_t* out_cursor = bitmap->buffer.data + bitmap->size_bits / 8;

  // First byte
  if ((out_i_cursor % 8) != 0) {
    int64_t n_partial_bits = _ArrowRoundUpToMultipleOf8(out_i_cursor) - out_i_cursor;
    for (int i = 0; i < n_partial_bits; i++) {
      ArrowBitSetTo(bitmap->buffer.data, out_i_cursor++, values[i]);
    }

    out_cursor++;
    values_cursor += n_partial_bits;
    n_remaining -= n_partial_bits;
  }

  // Middle bytes
  int64_t n_full_bytes = n_remaining / 8;
  for (int64_t i = 0; i < n_full_bytes; i++) {
    _ArrowBitmapPackInt32(values_cursor, out_cursor);
    values_cursor += 8;
    out_cursor++;
  }

  // Last byte
  out_i_cursor += n_full_bytes * 8;
  n_remaining -= n_full_bytes * 8;
  if (n_remaining > 0) {
    // Zero out the last byte
    *out_cursor = 0x00;
    for (int i = 0; i < n_remaining; i++) {
      ArrowBitSetTo(bitmap->buffer.data, out_i_cursor++, values_cursor[i]);
    }
    out_cursor++;
  }

  bitmap->size_bits += n_values;
  bitmap->buffer.size_bytes = out_cursor - bitmap->buffer.data;
}

static inline void ArrowBitmapReset(struct ArrowBitmap* bitmap) {
  ArrowBufferReset(&bitmap->buffer);
  bitmap->size_bits = 0;
}

#ifdef __cplusplus
}
#endif

#endif
// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#ifndef NANOARROW_ARRAY_INLINE_H_INCLUDED
#define NANOARROW_ARRAY_INLINE_H_INCLUDED

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>




#ifdef __cplusplus
extern "C" {
#endif

static inline struct ArrowBitmap* ArrowArrayValidityBitmap(struct ArrowArray* array) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;
  return &private_data->bitmap;
}

static inline struct ArrowBuffer* ArrowArrayBuffer(struct ArrowArray* array, int64_t i) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;
  switch (i) {
    case 0:
      return &private_data->bitmap.buffer;
    default:
      return private_data->buffers + i - 1;
  }
}

// We don't currently support the case of unions where type_id != child_index;
// however, these functions are used to keep track of where that assumption
// is made.
static inline int8_t _ArrowArrayUnionChildIndex(struct ArrowArray* array,
                                                int8_t type_id) {
  return type_id;
}

static inline int8_t _ArrowArrayUnionTypeId(struct ArrowArray* array,
                                            int8_t child_index) {
  return child_index;
}

static inline int8_t _ArrowParseUnionTypeIds(const char* type_ids, int8_t* out) {
  if (*type_ids == '\0') {
    return 0;
  }

  int32_t i = 0;
  long type_id;
  char* end_ptr;
  do {
    type_id = strtol(type_ids, &end_ptr, 10);
    if (end_ptr == type_ids || type_id < 0 || type_id > 127) {
      return -1;
    }

    if (out != NULL) {
      out[i] = (int8_t)type_id;
    }

    i++;

    type_ids = end_ptr;
    if (*type_ids == '\0') {
      return i;
    } else if (*type_ids != ',') {
      return -1;
    } else {
      type_ids++;
    }
  } while (1);

  return -1;
}

static inline int8_t _ArrowParsedUnionTypeIdsWillEqualChildIndices(const int8_t* type_ids,
                                                                   int64_t n_type_ids,
                                                                   int64_t n_children) {
  if (n_type_ids != n_children) {
    return 0;
  }

  for (int8_t i = 0; i < n_type_ids; i++) {
    if (type_ids[i] != i) {
      return 0;
    }
  }

  return 1;
}

static inline int8_t _ArrowUnionTypeIdsWillEqualChildIndices(const char* type_id_str,
                                                             int64_t n_children) {
  int8_t type_ids[128];
  int8_t n_type_ids = _ArrowParseUnionTypeIds(type_id_str, type_ids);
  return _ArrowParsedUnionTypeIdsWillEqualChildIndices(type_ids, n_type_ids, n_children);
}

static inline ArrowErrorCode ArrowArrayStartAppending(struct ArrowArray* array) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  switch (private_data->storage_type) {
    case NANOARROW_TYPE_UNINITIALIZED:
      return EINVAL;
    case NANOARROW_TYPE_SPARSE_UNION:
    case NANOARROW_TYPE_DENSE_UNION:
      // Note that this value could be -1 if the type_ids string was invalid
      if (private_data->union_type_id_is_child_index != 1) {
        return EINVAL;
      } else {
        break;
      }
    default:
      break;
  }
  if (private_data->storage_type == NANOARROW_TYPE_UNINITIALIZED) {
    return EINVAL;
  }

  // Initialize any data offset buffer with a single zero
  for (int i = 0; i < 3; i++) {
    if (private_data->layout.buffer_type[i] == NANOARROW_BUFFER_TYPE_DATA_OFFSET &&
        private_data->layout.element_size_bits[i] == 64) {
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt64(ArrowArrayBuffer(array, i), 0));
    } else if (private_data->layout.buffer_type[i] == NANOARROW_BUFFER_TYPE_DATA_OFFSET &&
               private_data->layout.element_size_bits[i] == 32) {
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt32(ArrowArrayBuffer(array, i), 0));
    }
  }

  // Start building any child arrays or dictionaries
  for (int64_t i = 0; i < array->n_children; i++) {
    NANOARROW_RETURN_NOT_OK(ArrowArrayStartAppending(array->children[i]));
  }

  if (array->dictionary != NULL) {
    NANOARROW_RETURN_NOT_OK(ArrowArrayStartAppending(array->dictionary));
  }

  return NANOARROW_OK;
}

static inline ArrowErrorCode ArrowArrayShrinkToFit(struct ArrowArray* array) {
  for (int64_t i = 0; i < 3; i++) {
    struct ArrowBuffer* buffer = ArrowArrayBuffer(array, i);
    NANOARROW_RETURN_NOT_OK(ArrowBufferResize(buffer, buffer->size_bytes, 1));
  }

  for (int64_t i = 0; i < array->n_children; i++) {
    NANOARROW_RETURN_NOT_OK(ArrowArrayShrinkToFit(array->children[i]));
  }

  if (array->dictionary != NULL) {
    NANOARROW_RETURN_NOT_OK(ArrowArrayShrinkToFit(array->dictionary));
  }

  return NANOARROW_OK;
}

static inline ArrowErrorCode _ArrowArrayAppendBits(struct ArrowArray* array,
                                                   int64_t buffer_i, uint8_t value,
                                                   int64_t n) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;
  struct ArrowBuffer* buffer = ArrowArrayBuffer(array, buffer_i);
  int64_t bytes_required =
      _ArrowRoundUpToMultipleOf8(private_data->layout.element_size_bits[buffer_i] *
                                 (array->length + 1)) /
      8;
  if (bytes_required > buffer->size_bytes) {
    NANOARROW_RETURN_NOT_OK(
        ArrowBufferAppendFill(buffer, 0, bytes_required - buffer->size_bytes));
  }

  ArrowBitsSetTo(buffer->data, array->length, n, value);
  return NANOARROW_OK;
}

static inline ArrowErrorCode _ArrowArrayAppendEmptyInternal(struct ArrowArray* array,
                                                            int64_t n, uint8_t is_valid) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  if (n == 0) {
    return NANOARROW_OK;
  }

  // Some type-specific handling
  switch (private_data->storage_type) {
    case NANOARROW_TYPE_NA:
      // (An empty value for a null array *is* a null)
      array->null_count += n;
      array->length += n;
      return NANOARROW_OK;

    case NANOARROW_TYPE_DENSE_UNION: {
      // Add one null to the first child and append n references to that child
      int8_t type_id = _ArrowArrayUnionTypeId(array, 0);
      NANOARROW_RETURN_NOT_OK(
          _ArrowArrayAppendEmptyInternal(array->children[0], 1, is_valid));
      NANOARROW_RETURN_NOT_OK(
          ArrowBufferAppendFill(ArrowArrayBuffer(array, 0), type_id, n));
      for (int64_t i = 0; i < n; i++) {
        NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt32(
            ArrowArrayBuffer(array, 1), (int32_t)array->children[0]->length - 1));
      }
      // For the purposes of array->null_count, union elements are never considered "null"
      // even if some children contain nulls.
      array->length += n;
      return NANOARROW_OK;
    }

    case NANOARROW_TYPE_SPARSE_UNION: {
      // Add n nulls to the first child and append n references to that child
      int8_t type_id = _ArrowArrayUnionTypeId(array, 0);
      NANOARROW_RETURN_NOT_OK(
          _ArrowArrayAppendEmptyInternal(array->children[0], n, is_valid));
      for (int64_t i = 1; i < array->n_children; i++) {
        NANOARROW_RETURN_NOT_OK(ArrowArrayAppendEmpty(array->children[i], n));
      }

      NANOARROW_RETURN_NOT_OK(
          ArrowBufferAppendFill(ArrowArrayBuffer(array, 0), type_id, n));
      // For the purposes of array->null_count, union elements are never considered "null"
      // even if some children contain nulls.
      array->length += n;
      return NANOARROW_OK;
    }

    case NANOARROW_TYPE_FIXED_SIZE_LIST:
      NANOARROW_RETURN_NOT_OK(ArrowArrayAppendEmpty(
          array->children[0], n * private_data->layout.child_size_elements));
      break;
    case NANOARROW_TYPE_STRUCT:
      for (int64_t i = 0; i < array->n_children; i++) {
        NANOARROW_RETURN_NOT_OK(ArrowArrayAppendEmpty(array->children[i], n));
      }
      break;

    default:
      break;
  }

  // Append n is_valid bits to the validity bitmap. If we haven't allocated a bitmap yet
  // and we need to append nulls, do it now.
  if (!is_valid && private_data->bitmap.buffer.data == NULL) {
    NANOARROW_RETURN_NOT_OK(ArrowBitmapReserve(&private_data->bitmap, array->length + n));
    ArrowBitmapAppendUnsafe(&private_data->bitmap, 1, array->length);
    ArrowBitmapAppendUnsafe(&private_data->bitmap, is_valid, n);
  } else if (private_data->bitmap.buffer.data != NULL) {
    NANOARROW_RETURN_NOT_OK(ArrowBitmapReserve(&private_data->bitmap, n));
    ArrowBitmapAppendUnsafe(&private_data->bitmap, is_valid, n);
  }

  // Add appropriate buffer fill
  struct ArrowBuffer* buffer;
  int64_t size_bytes;

  for (int i = 0; i < 3; i++) {
    buffer = ArrowArrayBuffer(array, i);
    size_bytes = private_data->layout.element_size_bits[i] / 8;

    switch (private_data->layout.buffer_type[i]) {
      case NANOARROW_BUFFER_TYPE_NONE:
      case NANOARROW_BUFFER_TYPE_VALIDITY:
        continue;
      case NANOARROW_BUFFER_TYPE_DATA_OFFSET:
        // Append the current value at the end of the offset buffer for each element
        NANOARROW_RETURN_NOT_OK(ArrowBufferReserve(buffer, size_bytes * n));

        for (int64_t j = 0; j < n; j++) {
          ArrowBufferAppendUnsafe(buffer, buffer->data + size_bytes * (array->length + j),
                                  size_bytes);
        }

        // Skip the data buffer
        i++;
        continue;
      case NANOARROW_BUFFER_TYPE_DATA:
        // Zero out the next bit of memory
        if (private_data->layout.element_size_bits[i] % 8 == 0) {
          NANOARROW_RETURN_NOT_OK(ArrowBufferAppendFill(buffer, 0, size_bytes * n));
        } else {
          NANOARROW_RETURN_NOT_OK(_ArrowArrayAppendBits(array, i, 0, n));
        }
        continue;

      case NANOARROW_BUFFER_TYPE_TYPE_ID:
      case NANOARROW_BUFFER_TYPE_UNION_OFFSET:
        // These cases return above
        return EINVAL;
    }
  }

  array->length += n;
  array->null_count += n * !is_valid;
  return NANOARROW_OK;
}

static inline ArrowErrorCode ArrowArrayAppendNull(struct ArrowArray* array, int64_t n) {
  return _ArrowArrayAppendEmptyInternal(array, n, 0);
}

static inline ArrowErrorCode ArrowArrayAppendEmpty(struct ArrowArray* array, int64_t n) {
  return _ArrowArrayAppendEmptyInternal(array, n, 1);
}

static inline ArrowErrorCode ArrowArrayAppendInt(struct ArrowArray* array,
                                                 int64_t value) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  struct ArrowBuffer* data_buffer = ArrowArrayBuffer(array, 1);

  switch (private_data->storage_type) {
    case NANOARROW_TYPE_INT64:
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppend(data_buffer, &value, sizeof(int64_t)));
      break;
    case NANOARROW_TYPE_INT32:
      _NANOARROW_CHECK_RANGE(value, INT32_MIN, INT32_MAX);
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt32(data_buffer, (int32_t)value));
      break;
    case NANOARROW_TYPE_INT16:
      _NANOARROW_CHECK_RANGE(value, INT16_MIN, INT16_MAX);
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt16(data_buffer, (int16_t)value));
      break;
    case NANOARROW_TYPE_INT8:
      _NANOARROW_CHECK_RANGE(value, INT8_MIN, INT8_MAX);
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt8(data_buffer, (int8_t)value));
      break;
    case NANOARROW_TYPE_UINT64:
    case NANOARROW_TYPE_UINT32:
    case NANOARROW_TYPE_UINT16:
    case NANOARROW_TYPE_UINT8:
      _NANOARROW_CHECK_RANGE(value, 0, INT64_MAX);
      return ArrowArrayAppendUInt(array, value);
    case NANOARROW_TYPE_DOUBLE:
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendDouble(data_buffer, (double)value));
      break;
    case NANOARROW_TYPE_FLOAT:
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendFloat(data_buffer, (float)value));
      break;
    case NANOARROW_TYPE_BOOL:
      NANOARROW_RETURN_NOT_OK(_ArrowArrayAppendBits(array, 1, value != 0, 1));
      break;
    default:
      return EINVAL;
  }

  if (private_data->bitmap.buffer.data != NULL) {
    NANOARROW_RETURN_NOT_OK(ArrowBitmapAppend(ArrowArrayValidityBitmap(array), 1, 1));
  }

  array->length++;
  return NANOARROW_OK;
}

static inline ArrowErrorCode ArrowArrayAppendUInt(struct ArrowArray* array,
                                                  uint64_t value) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  struct ArrowBuffer* data_buffer = ArrowArrayBuffer(array, 1);

  switch (private_data->storage_type) {
    case NANOARROW_TYPE_UINT64:
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppend(data_buffer, &value, sizeof(uint64_t)));
      break;
    case NANOARROW_TYPE_UINT32:
      _NANOARROW_CHECK_UPPER_LIMIT(value, UINT32_MAX);
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendUInt32(data_buffer, (uint32_t)value));
      break;
    case NANOARROW_TYPE_UINT16:
      _NANOARROW_CHECK_UPPER_LIMIT(value, UINT16_MAX);
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendUInt16(data_buffer, (uint16_t)value));
      break;
    case NANOARROW_TYPE_UINT8:
      _NANOARROW_CHECK_UPPER_LIMIT(value, UINT8_MAX);
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendUInt8(data_buffer, (uint8_t)value));
      break;
    case NANOARROW_TYPE_INT64:
    case NANOARROW_TYPE_INT32:
    case NANOARROW_TYPE_INT16:
    case NANOARROW_TYPE_INT8:
      _NANOARROW_CHECK_UPPER_LIMIT(value, INT64_MAX);
      return ArrowArrayAppendInt(array, value);
    case NANOARROW_TYPE_DOUBLE:
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendDouble(data_buffer, (double)value));
      break;
    case NANOARROW_TYPE_FLOAT:
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendFloat(data_buffer, (float)value));
      break;
    case NANOARROW_TYPE_BOOL:
      NANOARROW_RETURN_NOT_OK(_ArrowArrayAppendBits(array, 1, value != 0, 1));
      break;
    default:
      return EINVAL;
  }

  if (private_data->bitmap.buffer.data != NULL) {
    NANOARROW_RETURN_NOT_OK(ArrowBitmapAppend(ArrowArrayValidityBitmap(array), 1, 1));
  }

  array->length++;
  return NANOARROW_OK;
}

static inline ArrowErrorCode ArrowArrayAppendDouble(struct ArrowArray* array,
                                                    double value) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  struct ArrowBuffer* data_buffer = ArrowArrayBuffer(array, 1);

  switch (private_data->storage_type) {
    case NANOARROW_TYPE_DOUBLE:
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppend(data_buffer, &value, sizeof(double)));
      break;
    case NANOARROW_TYPE_FLOAT:
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendFloat(data_buffer, (float)value));
      break;
    default:
      return EINVAL;
  }

  if (private_data->bitmap.buffer.data != NULL) {
    NANOARROW_RETURN_NOT_OK(ArrowBitmapAppend(ArrowArrayValidityBitmap(array), 1, 1));
  }

  array->length++;
  return NANOARROW_OK;
}

static inline ArrowErrorCode ArrowArrayAppendBytes(struct ArrowArray* array,
                                                   struct ArrowBufferView value) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  struct ArrowBuffer* offset_buffer = ArrowArrayBuffer(array, 1);
  struct ArrowBuffer* data_buffer = ArrowArrayBuffer(
      array, 1 + (private_data->storage_type != NANOARROW_TYPE_FIXED_SIZE_BINARY));
  int32_t offset;
  int64_t large_offset;
  int64_t fixed_size_bytes = private_data->layout.element_size_bits[1] / 8;

  switch (private_data->storage_type) {
    case NANOARROW_TYPE_STRING:
    case NANOARROW_TYPE_BINARY:
      offset = ((int32_t*)offset_buffer->data)[array->length];
      if ((offset + value.size_bytes) > INT32_MAX) {
        return EINVAL;
      }

      offset += (int32_t)value.size_bytes;
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppend(offset_buffer, &offset, sizeof(int32_t)));
      NANOARROW_RETURN_NOT_OK(
          ArrowBufferAppend(data_buffer, value.data.data, value.size_bytes));
      break;

    case NANOARROW_TYPE_LARGE_STRING:
    case NANOARROW_TYPE_LARGE_BINARY:
      large_offset = ((int64_t*)offset_buffer->data)[array->length];
      large_offset += value.size_bytes;
      NANOARROW_RETURN_NOT_OK(
          ArrowBufferAppend(offset_buffer, &large_offset, sizeof(int64_t)));
      NANOARROW_RETURN_NOT_OK(
          ArrowBufferAppend(data_buffer, value.data.data, value.size_bytes));
      break;

    case NANOARROW_TYPE_FIXED_SIZE_BINARY:
      if (value.size_bytes != fixed_size_bytes) {
        return EINVAL;
      }

      NANOARROW_RETURN_NOT_OK(
          ArrowBufferAppend(data_buffer, value.data.data, value.size_bytes));
      break;
    default:
      return EINVAL;
  }

  if (private_data->bitmap.buffer.data != NULL) {
    NANOARROW_RETURN_NOT_OK(ArrowBitmapAppend(ArrowArrayValidityBitmap(array), 1, 1));
  }

  array->length++;
  return NANOARROW_OK;
}

static inline ArrowErrorCode ArrowArrayAppendString(struct ArrowArray* array,
                                                    struct ArrowStringView value) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  struct ArrowBufferView buffer_view;
  buffer_view.data.data = value.data;
  buffer_view.size_bytes = value.size_bytes;

  switch (private_data->storage_type) {
    case NANOARROW_TYPE_STRING:
    case NANOARROW_TYPE_LARGE_STRING:
    case NANOARROW_TYPE_BINARY:
    case NANOARROW_TYPE_LARGE_BINARY:
      return ArrowArrayAppendBytes(array, buffer_view);
    default:
      return EINVAL;
  }
}

static inline ArrowErrorCode ArrowArrayAppendInterval(struct ArrowArray* array,
                                                      struct ArrowInterval* value) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  struct ArrowBuffer* data_buffer = ArrowArrayBuffer(array, 1);

  switch (private_data->storage_type) {
    case NANOARROW_TYPE_INTERVAL_MONTHS: {
      if (value->type != NANOARROW_TYPE_INTERVAL_MONTHS) {
        return EINVAL;
      }

      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt32(data_buffer, value->months));
      break;
    }
    case NANOARROW_TYPE_INTERVAL_DAY_TIME: {
      if (value->type != NANOARROW_TYPE_INTERVAL_DAY_TIME) {
        return EINVAL;
      }

      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt32(data_buffer, value->days));
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt32(data_buffer, value->ms));
      break;
    }
    case NANOARROW_TYPE_INTERVAL_MONTH_DAY_NANO: {
      if (value->type != NANOARROW_TYPE_INTERVAL_MONTH_DAY_NANO) {
        return EINVAL;
      }

      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt32(data_buffer, value->months));
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt32(data_buffer, value->days));
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt64(data_buffer, value->ns));
      break;
    }
    default:
      return EINVAL;
  }

  array->length++;
  return NANOARROW_OK;
}

static inline ArrowErrorCode ArrowArrayAppendDecimal(struct ArrowArray* array,
                                                     struct ArrowDecimal* value) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;
  struct ArrowBuffer* data_buffer = ArrowArrayBuffer(array, 1);

  switch (private_data->storage_type) {
    case NANOARROW_TYPE_DECIMAL128:
      if (value->n_words != 2) {
        return EINVAL;
      } else {
        NANOARROW_RETURN_NOT_OK(
            ArrowBufferAppend(data_buffer, value->words, 2 * sizeof(uint64_t)));
        break;
      }
    case NANOARROW_TYPE_DECIMAL256:
      if (value->n_words != 4) {
        return EINVAL;
      } else {
        NANOARROW_RETURN_NOT_OK(
            ArrowBufferAppend(data_buffer, value->words, 4 * sizeof(uint64_t)));
        break;
      }
    default:
      return EINVAL;
  }

  if (private_data->bitmap.buffer.data != NULL) {
    NANOARROW_RETURN_NOT_OK(ArrowBitmapAppend(ArrowArrayValidityBitmap(array), 1, 1));
  }

  array->length++;
  return NANOARROW_OK;
}

static inline ArrowErrorCode ArrowArrayFinishElement(struct ArrowArray* array) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  int64_t child_length;

  switch (private_data->storage_type) {
    case NANOARROW_TYPE_LIST:
    case NANOARROW_TYPE_MAP:
      child_length = array->children[0]->length;
      if (child_length > INT32_MAX) {
        return EINVAL;
      }
      NANOARROW_RETURN_NOT_OK(
          ArrowBufferAppendInt32(ArrowArrayBuffer(array, 1), (int32_t)child_length));
      break;
    case NANOARROW_TYPE_LARGE_LIST:
      child_length = array->children[0]->length;
      NANOARROW_RETURN_NOT_OK(
          ArrowBufferAppendInt64(ArrowArrayBuffer(array, 1), child_length));
      break;
    case NANOARROW_TYPE_FIXED_SIZE_LIST:
      child_length = array->children[0]->length;
      if (child_length !=
          ((array->length + 1) * private_data->layout.child_size_elements)) {
        return EINVAL;
      }
      break;
    case NANOARROW_TYPE_STRUCT:
      for (int64_t i = 0; i < array->n_children; i++) {
        child_length = array->children[i]->length;
        if (child_length != (array->length + 1)) {
          return EINVAL;
        }
      }
      break;
    default:
      return EINVAL;
  }

  if (private_data->bitmap.buffer.data != NULL) {
    NANOARROW_RETURN_NOT_OK(ArrowBitmapAppend(ArrowArrayValidityBitmap(array), 1, 1));
  }

  array->length++;
  return NANOARROW_OK;
}

static inline ArrowErrorCode ArrowArrayFinishUnionElement(struct ArrowArray* array,
                                                          int8_t type_id) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  int64_t child_index = _ArrowArrayUnionChildIndex(array, type_id);
  if (child_index < 0 || child_index >= array->n_children) {
    return EINVAL;
  }

  switch (private_data->storage_type) {
    case NANOARROW_TYPE_DENSE_UNION:
      // Append the target child length to the union offsets buffer
      _NANOARROW_CHECK_RANGE(array->children[child_index]->length, 0, INT32_MAX);
      NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt32(
          ArrowArrayBuffer(array, 1), (int32_t)array->children[child_index]->length - 1));
      break;
    case NANOARROW_TYPE_SPARSE_UNION:
      // Append one empty to any non-target column that isn't already the right length
      // or abort if appending a null will result in a column with invalid length
      for (int64_t i = 0; i < array->n_children; i++) {
        if (i == child_index || array->children[i]->length == (array->length + 1)) {
          continue;
        }

        if (array->children[i]->length != array->length) {
          return EINVAL;
        }

        NANOARROW_RETURN_NOT_OK(ArrowArrayAppendEmpty(array->children[i], 1));
      }

      break;
    default:
      return EINVAL;
  }

  // Write to the type_ids buffer
  NANOARROW_RETURN_NOT_OK(
      ArrowBufferAppendInt8(ArrowArrayBuffer(array, 0), (int8_t)type_id));
  array->length++;
  return NANOARROW_OK;
}

static inline void ArrowArrayViewMove(struct ArrowArrayView* src,
                                      struct ArrowArrayView* dst) {
  memcpy(dst, src, sizeof(struct ArrowArrayView));
  ArrowArrayViewInitFromType(src, NANOARROW_TYPE_UNINITIALIZED);
}

static inline int8_t ArrowArrayViewIsNull(struct ArrowArrayView* array_view, int64_t i) {
  const uint8_t* validity_buffer = array_view->buffer_views[0].data.as_uint8;
  i += array_view->offset;
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_NA:
      return 0x01;
    case NANOARROW_TYPE_DENSE_UNION:
    case NANOARROW_TYPE_SPARSE_UNION:
      // Unions are "never null" in Arrow land
      return 0x00;
    default:
      return validity_buffer != NULL && !ArrowBitGet(validity_buffer, i);
  }
}

static inline int8_t ArrowArrayViewUnionTypeId(struct ArrowArrayView* array_view,
                                               int64_t i) {
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_DENSE_UNION:
    case NANOARROW_TYPE_SPARSE_UNION:
      return array_view->buffer_views[0].data.as_int8[i];
    default:
      return -1;
  }
}

static inline int8_t ArrowArrayViewUnionChildIndex(struct ArrowArrayView* array_view,
                                                   int64_t i) {
  int8_t type_id = ArrowArrayViewUnionTypeId(array_view, i);
  if (array_view->union_type_id_map == NULL) {
    return type_id;
  } else {
    return array_view->union_type_id_map[type_id];
  }
}

static inline int64_t ArrowArrayViewUnionChildOffset(struct ArrowArrayView* array_view,
                                                     int64_t i) {
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_DENSE_UNION:
      return array_view->buffer_views[1].data.as_int32[i];
    case NANOARROW_TYPE_SPARSE_UNION:
      return i;
    default:
      return -1;
  }
}

static inline int64_t ArrowArrayViewListChildOffset(struct ArrowArrayView* array_view,
                                                    int64_t i) {
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_LIST:
      return array_view->buffer_views[1].data.as_int32[i];
    case NANOARROW_TYPE_LARGE_LIST:
      return array_view->buffer_views[1].data.as_int64[i];
    default:
      return -1;
  }
}

static inline int64_t ArrowArrayViewGetIntUnsafe(struct ArrowArrayView* array_view,
                                                 int64_t i) {
  struct ArrowBufferView* data_view = &array_view->buffer_views[1];
  i += array_view->offset;
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_INT64:
      return data_view->data.as_int64[i];
    case NANOARROW_TYPE_UINT64:
      return data_view->data.as_uint64[i];
    case NANOARROW_TYPE_INT32:
      return data_view->data.as_int32[i];
    case NANOARROW_TYPE_UINT32:
      return data_view->data.as_uint32[i];
    case NANOARROW_TYPE_INT16:
      return data_view->data.as_int16[i];
    case NANOARROW_TYPE_UINT16:
      return data_view->data.as_uint16[i];
    case NANOARROW_TYPE_INT8:
      return data_view->data.as_int8[i];
    case NANOARROW_TYPE_UINT8:
      return data_view->data.as_uint8[i];
    case NANOARROW_TYPE_DOUBLE:
      return (int64_t)data_view->data.as_double[i];
    case NANOARROW_TYPE_FLOAT:
      return (int64_t)data_view->data.as_float[i];
    case NANOARROW_TYPE_BOOL:
      return ArrowBitGet(data_view->data.as_uint8, i);
    default:
      return INT64_MAX;
  }
}

static inline uint64_t ArrowArrayViewGetUIntUnsafe(struct ArrowArrayView* array_view,
                                                   int64_t i) {
  i += array_view->offset;
  struct ArrowBufferView* data_view = &array_view->buffer_views[1];
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_INT64:
      return data_view->data.as_int64[i];
    case NANOARROW_TYPE_UINT64:
      return data_view->data.as_uint64[i];
    case NANOARROW_TYPE_INT32:
      return data_view->data.as_int32[i];
    case NANOARROW_TYPE_UINT32:
      return data_view->data.as_uint32[i];
    case NANOARROW_TYPE_INT16:
      return data_view->data.as_int16[i];
    case NANOARROW_TYPE_UINT16:
      return data_view->data.as_uint16[i];
    case NANOARROW_TYPE_INT8:
      return data_view->data.as_int8[i];
    case NANOARROW_TYPE_UINT8:
      return data_view->data.as_uint8[i];
    case NANOARROW_TYPE_DOUBLE:
      return (uint64_t)data_view->data.as_double[i];
    case NANOARROW_TYPE_FLOAT:
      return (uint64_t)data_view->data.as_float[i];
    case NANOARROW_TYPE_BOOL:
      return ArrowBitGet(data_view->data.as_uint8, i);
    default:
      return UINT64_MAX;
  }
}

static inline double ArrowArrayViewGetDoubleUnsafe(struct ArrowArrayView* array_view,
                                                   int64_t i) {
  i += array_view->offset;
  struct ArrowBufferView* data_view = &array_view->buffer_views[1];
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_INT64:
      return (double)data_view->data.as_int64[i];
    case NANOARROW_TYPE_UINT64:
      return (double)data_view->data.as_uint64[i];
    case NANOARROW_TYPE_INT32:
      return data_view->data.as_int32[i];
    case NANOARROW_TYPE_UINT32:
      return data_view->data.as_uint32[i];
    case NANOARROW_TYPE_INT16:
      return data_view->data.as_int16[i];
    case NANOARROW_TYPE_UINT16:
      return data_view->data.as_uint16[i];
    case NANOARROW_TYPE_INT8:
      return data_view->data.as_int8[i];
    case NANOARROW_TYPE_UINT8:
      return data_view->data.as_uint8[i];
    case NANOARROW_TYPE_DOUBLE:
      return data_view->data.as_double[i];
    case NANOARROW_TYPE_FLOAT:
      return data_view->data.as_float[i];
    case NANOARROW_TYPE_BOOL:
      return ArrowBitGet(data_view->data.as_uint8, i);
    default:
      return DBL_MAX;
  }
}

static inline struct ArrowStringView ArrowArrayViewGetStringUnsafe(
    struct ArrowArrayView* array_view, int64_t i) {
  i += array_view->offset;
  struct ArrowBufferView* offsets_view = &array_view->buffer_views[1];
  const char* data_view = array_view->buffer_views[2].data.as_char;

  struct ArrowStringView view;
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_STRING:
    case NANOARROW_TYPE_BINARY:
      view.data = data_view + offsets_view->data.as_int32[i];
      view.size_bytes =
          offsets_view->data.as_int32[i + 1] - offsets_view->data.as_int32[i];
      break;
    case NANOARROW_TYPE_LARGE_STRING:
    case NANOARROW_TYPE_LARGE_BINARY:
      view.data = data_view + offsets_view->data.as_int64[i];
      view.size_bytes =
          offsets_view->data.as_int64[i + 1] - offsets_view->data.as_int64[i];
      break;
    case NANOARROW_TYPE_FIXED_SIZE_BINARY:
      view.size_bytes = array_view->layout.element_size_bits[1] / 8;
      view.data = array_view->buffer_views[1].data.as_char + (i * view.size_bytes);
      break;
    default:
      view.data = NULL;
      view.size_bytes = 0;
      break;
  }

  return view;
}

static inline struct ArrowBufferView ArrowArrayViewGetBytesUnsafe(
    struct ArrowArrayView* array_view, int64_t i) {
  i += array_view->offset;
  struct ArrowBufferView* offsets_view = &array_view->buffer_views[1];
  const uint8_t* data_view = array_view->buffer_views[2].data.as_uint8;

  struct ArrowBufferView view;
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_STRING:
    case NANOARROW_TYPE_BINARY:
      view.size_bytes =
          offsets_view->data.as_int32[i + 1] - offsets_view->data.as_int32[i];
      view.data.as_uint8 = data_view + offsets_view->data.as_int32[i];
      break;
    case NANOARROW_TYPE_LARGE_STRING:
    case NANOARROW_TYPE_LARGE_BINARY:
      view.size_bytes =
          offsets_view->data.as_int64[i + 1] - offsets_view->data.as_int64[i];
      view.data.as_uint8 = data_view + offsets_view->data.as_int64[i];
      break;
    case NANOARROW_TYPE_FIXED_SIZE_BINARY:
      view.size_bytes = array_view->layout.element_size_bits[1] / 8;
      view.data.as_uint8 =
          array_view->buffer_views[1].data.as_uint8 + (i * view.size_bytes);
      break;
    default:
      view.data.data = NULL;
      view.size_bytes = 0;
      break;
  }

  return view;
}

static inline void ArrowArrayViewGetIntervalUnsafe(struct ArrowArrayView* array_view,
                                                   int64_t i, struct ArrowInterval* out) {
  const uint8_t* data_view = array_view->buffer_views[1].data.as_uint8;
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_INTERVAL_MONTHS: {
      const size_t size = sizeof(int32_t);
      memcpy(&out->months, data_view + i * size, sizeof(int32_t));
      break;
    }
    case NANOARROW_TYPE_INTERVAL_DAY_TIME: {
      const size_t size = sizeof(int32_t) + sizeof(int32_t);
      memcpy(&out->days, data_view + i * size, sizeof(int32_t));
      memcpy(&out->ms, data_view + i * size + 4, sizeof(int32_t));
      break;
    }
    case NANOARROW_TYPE_INTERVAL_MONTH_DAY_NANO: {
      const size_t size = sizeof(int32_t) + sizeof(int32_t) + sizeof(int64_t);
      memcpy(&out->months, data_view + i * size, sizeof(int32_t));
      memcpy(&out->days, data_view + i * size + 4, sizeof(int32_t));
      memcpy(&out->ns, data_view + i * size + 8, sizeof(int64_t));
      break;
    }
    default:
      break;
  }
}

static inline void ArrowArrayViewGetDecimalUnsafe(struct ArrowArrayView* array_view,
                                                  int64_t i, struct ArrowDecimal* out) {
  i += array_view->offset;
  const uint8_t* data_view = array_view->buffer_views[1].data.as_uint8;
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_DECIMAL128:
      ArrowDecimalSetBytes(out, data_view + (i * 16));
      break;
    case NANOARROW_TYPE_DECIMAL256:
      ArrowDecimalSetBytes(out, data_view + (i * 32));
      break;
    default:
      memset(out->words, 0, sizeof(out->words));
      break;
  }
}

#ifdef __cplusplus
}
#endif

#endif
