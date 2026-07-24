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

#pragma once

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arrow/c/abi.h"

#define ARROW_C_ASSERT(condition, msg)                          \
  do {                                                          \
    if (!(condition)) {                                         \
      fprintf(stderr, "%s:%d:: %s", __FILE__, __LINE__, (msg)); \
      abort();                                                  \
    }                                                           \
  } while (0)

#ifdef __cplusplus
extern "C" {
#endif

/// Query whether the C schema is released
inline int ArrowSchemaIsReleased(const struct ArrowSchema* schema) {
  return schema->release == NULL;
}

/// Mark the C schema released (for use in release callbacks)
inline void ArrowSchemaMarkReleased(struct ArrowSchema* schema) {
  schema->release = NULL;
}

/// Move the C schema from `src` to `dest`
///
/// Note `dest` must *not* point to a valid schema already, otherwise there
/// will be a memory leak.
inline void ArrowSchemaMove(struct ArrowSchema* src, struct ArrowSchema* dest) {
  assert(dest != src);
  assert(!ArrowSchemaIsReleased(src));
  memcpy(dest, src, sizeof(struct ArrowSchema));
  ArrowSchemaMarkReleased(src);
}

/// Release the C schema, if necessary, by calling its release callback
inline void ArrowSchemaRelease(struct ArrowSchema* schema) {
  if (!ArrowSchemaIsReleased(schema)) {
    schema->release(schema);
    ARROW_C_ASSERT(ArrowSchemaIsReleased(schema),
                   "ArrowSchemaRelease did not cleanup release callback");
  }
}

/// Query whether the C array is released
inline int ArrowArrayIsReleased(const struct ArrowArray* array) {
  return array->release == NULL;
}

inline int ArrowDeviceArrayIsReleased(const struct ArrowDeviceArray* array) {
  return ArrowArrayIsReleased(&array->array);
}

/// Mark the C array released (for use in release callbacks)
inline void ArrowArrayMarkReleased(struct ArrowArray* array) { array->release = NULL; }

inline void ArrowDeviceArrayMarkReleased(struct ArrowDeviceArray* array) {
  ArrowArrayMarkReleased(&array->array);
}

/// Move the C array from `src` to `dest`
///
/// Note `dest` must *not* point to a valid array already, otherwise there
/// will be a memory leak.
inline void ArrowArrayMove(struct ArrowArray* src, struct ArrowArray* dest) {
  assert(dest != src);
  assert(!ArrowArrayIsReleased(src));
  memcpy(dest, src, sizeof(struct ArrowArray));
  ArrowArrayMarkReleased(src);
}

inline void ArrowDeviceArrayMove(struct ArrowDeviceArray* src,
                                 struct ArrowDeviceArray* dest) {
  assert(dest != src);
  assert(!ArrowDeviceArrayIsReleased(src));
  memcpy(dest, src, sizeof(struct ArrowDeviceArray));
  ArrowDeviceArrayMarkReleased(src);
}

/// Release the C array, if necessary, by calling its release callback
inline void ArrowArrayRelease(struct ArrowArray* array) {
  if (!ArrowArrayIsReleased(array)) {
    array->release(array);
    ARROW_C_ASSERT(ArrowArrayIsReleased(array),
                   "ArrowArrayRelease did not cleanup release callback");
  }
}

inline void ArrowDeviceArrayRelease(struct ArrowDeviceArray* array) {
  if (!ArrowDeviceArrayIsReleased(array)) {
    array->array.release(&array->array);
    ARROW_C_ASSERT(ArrowDeviceArrayIsReleased(array),
                   "ArrowDeviceArrayRelease did not cleanup release callback");
  }
}

/// Query whether the C array stream is released
inline int ArrowArrayStreamIsReleased(const struct ArrowArrayStream* stream) {
  return stream->release == NULL;
}

inline int ArrowDeviceArrayStreamIsReleased(const struct ArrowDeviceArrayStream* stream) {
  return stream->release == NULL;
}

/// Mark the C array stream released (for use in release callbacks)
inline void ArrowArrayStreamMarkReleased(struct ArrowArrayStream* stream) {
  stream->release = NULL;
}

inline void ArrowDeviceArrayStreamMarkReleased(struct ArrowDeviceArrayStream* stream) {
  stream->release = NULL;
}

/// Move the C array stream from `src` to `dest`
///
/// Note `dest` must *not* point to a valid stream already, otherwise there
/// will be a memory leak.
inline void ArrowArrayStreamMove(struct ArrowArrayStream* src,
                                 struct ArrowArrayStream* dest) {
  assert(dest != src);
  assert(!ArrowArrayStreamIsReleased(src));
  memcpy(dest, src, sizeof(struct ArrowArrayStream));
  ArrowArrayStreamMarkReleased(src);
}

inline void ArrowDeviceArrayStreamMove(struct ArrowDeviceArrayStream* src,
                                       struct ArrowDeviceArrayStream* dest) {
  assert(dest != src);
  assert(!ArrowDeviceArrayStreamIsReleased(src));
  memcpy(dest, src, sizeof(struct ArrowDeviceArrayStream));
  ArrowDeviceArrayStreamMarkReleased(src);
}

/// Release the C array stream, if necessary, by calling its release callback
inline void ArrowArrayStreamRelease(struct ArrowArrayStream* stream) {
  if (!ArrowArrayStreamIsReleased(stream)) {
    stream->release(stream);
    ARROW_C_ASSERT(ArrowArrayStreamIsReleased(stream),
                   "ArrowArrayStreamRelease did not cleanup release callback");
  }
}

inline void ArrowDeviceArrayStreamRelease(struct ArrowDeviceArrayStream* stream) {
  if (!ArrowDeviceArrayStreamIsReleased(stream)) {
    stream->release(stream);
    ARROW_C_ASSERT(ArrowDeviceArrayStreamIsReleased(stream),
                   "ArrowDeviceArrayStreamRelease did not cleanup release callback");
  }
}

#ifdef __cplusplus
}
#endif
