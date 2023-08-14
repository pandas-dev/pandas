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

#include <errno.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nanoarrow.h"

const char* ArrowNanoarrowVersion(void) { return NANOARROW_VERSION; }

int ArrowNanoarrowVersionInt(void) { return NANOARROW_VERSION_INT; }

int ArrowErrorSet(struct ArrowError* error, const char* fmt, ...) {
  if (error == NULL) {
    return NANOARROW_OK;
  }

  memset(error->message, 0, sizeof(error->message));

  va_list args;
  va_start(args, fmt);
  int chars_needed = vsnprintf(error->message, sizeof(error->message), fmt, args);
  va_end(args);

  if (chars_needed < 0) {
    return EINVAL;
  } else if (((size_t)chars_needed) >= sizeof(error->message)) {
    return ERANGE;
  } else {
    return NANOARROW_OK;
  }
}

const char* ArrowErrorMessage(struct ArrowError* error) {
  if (error == NULL) {
    return "";
  } else {
    return error->message;
  }
}

void ArrowLayoutInit(struct ArrowLayout* layout, enum ArrowType storage_type) {
  layout->buffer_type[0] = NANOARROW_BUFFER_TYPE_VALIDITY;
  layout->buffer_data_type[0] = NANOARROW_TYPE_BOOL;
  layout->buffer_type[1] = NANOARROW_BUFFER_TYPE_DATA;
  layout->buffer_data_type[1] = storage_type;
  layout->buffer_type[2] = NANOARROW_BUFFER_TYPE_NONE;
  layout->buffer_data_type[2] = NANOARROW_TYPE_UNINITIALIZED;

  layout->element_size_bits[0] = 1;
  layout->element_size_bits[1] = 0;
  layout->element_size_bits[2] = 0;

  layout->child_size_elements = 0;

  switch (storage_type) {
    case NANOARROW_TYPE_UNINITIALIZED:
    case NANOARROW_TYPE_NA:
      layout->buffer_type[0] = NANOARROW_BUFFER_TYPE_NONE;
      layout->buffer_data_type[0] = NANOARROW_TYPE_UNINITIALIZED;
      layout->buffer_type[1] = NANOARROW_BUFFER_TYPE_NONE;
      layout->buffer_data_type[1] = NANOARROW_TYPE_UNINITIALIZED;
      layout->element_size_bits[0] = 0;
      break;

    case NANOARROW_TYPE_LIST:
    case NANOARROW_TYPE_MAP:
      layout->buffer_type[1] = NANOARROW_BUFFER_TYPE_DATA_OFFSET;
      layout->buffer_data_type[1] = NANOARROW_TYPE_INT32;
      layout->element_size_bits[1] = 32;
      break;

    case NANOARROW_TYPE_LARGE_LIST:
      layout->buffer_type[1] = NANOARROW_BUFFER_TYPE_DATA_OFFSET;
      layout->buffer_data_type[1] = NANOARROW_TYPE_INT64;
      layout->element_size_bits[1] = 64;
      break;

    case NANOARROW_TYPE_STRUCT:
    case NANOARROW_TYPE_FIXED_SIZE_LIST:
      layout->buffer_type[1] = NANOARROW_BUFFER_TYPE_NONE;
      layout->buffer_data_type[1] = NANOARROW_TYPE_UNINITIALIZED;
      break;

    case NANOARROW_TYPE_BOOL:
      layout->element_size_bits[1] = 1;
      break;

    case NANOARROW_TYPE_UINT8:
    case NANOARROW_TYPE_INT8:
      layout->element_size_bits[1] = 8;
      break;

    case NANOARROW_TYPE_UINT16:
    case NANOARROW_TYPE_INT16:
    case NANOARROW_TYPE_HALF_FLOAT:
      layout->element_size_bits[1] = 16;
      break;

    case NANOARROW_TYPE_UINT32:
    case NANOARROW_TYPE_INT32:
    case NANOARROW_TYPE_FLOAT:
      layout->element_size_bits[1] = 32;
      break;
    case NANOARROW_TYPE_INTERVAL_MONTHS:
      layout->buffer_data_type[1] = NANOARROW_TYPE_INT32;
      layout->element_size_bits[1] = 32;
      break;

    case NANOARROW_TYPE_UINT64:
    case NANOARROW_TYPE_INT64:
    case NANOARROW_TYPE_DOUBLE:
    case NANOARROW_TYPE_INTERVAL_DAY_TIME:
      layout->element_size_bits[1] = 64;
      break;

    case NANOARROW_TYPE_DECIMAL128:
    case NANOARROW_TYPE_INTERVAL_MONTH_DAY_NANO:
      layout->element_size_bits[1] = 128;
      break;

    case NANOARROW_TYPE_DECIMAL256:
      layout->element_size_bits[1] = 256;
      break;

    case NANOARROW_TYPE_FIXED_SIZE_BINARY:
      layout->buffer_data_type[1] = NANOARROW_TYPE_BINARY;
      break;

    case NANOARROW_TYPE_DENSE_UNION:
      layout->buffer_type[0] = NANOARROW_BUFFER_TYPE_TYPE_ID;
      layout->buffer_data_type[0] = NANOARROW_TYPE_INT8;
      layout->element_size_bits[0] = 8;
      layout->buffer_type[1] = NANOARROW_BUFFER_TYPE_UNION_OFFSET;
      layout->buffer_data_type[1] = NANOARROW_TYPE_INT32;
      layout->element_size_bits[1] = 32;
      break;

    case NANOARROW_TYPE_SPARSE_UNION:
      layout->buffer_type[0] = NANOARROW_BUFFER_TYPE_TYPE_ID;
      layout->buffer_data_type[0] = NANOARROW_TYPE_INT8;
      layout->element_size_bits[0] = 8;
      layout->buffer_type[1] = NANOARROW_BUFFER_TYPE_NONE;
      layout->buffer_data_type[1] = NANOARROW_TYPE_UNINITIALIZED;
      break;

    case NANOARROW_TYPE_STRING:
    case NANOARROW_TYPE_BINARY:
      layout->buffer_type[1] = NANOARROW_BUFFER_TYPE_DATA_OFFSET;
      layout->buffer_data_type[1] = NANOARROW_TYPE_INT32;
      layout->element_size_bits[1] = 32;
      layout->buffer_type[2] = NANOARROW_BUFFER_TYPE_DATA;
      layout->buffer_data_type[2] = storage_type;
      break;

    case NANOARROW_TYPE_LARGE_STRING:
      layout->buffer_type[1] = NANOARROW_BUFFER_TYPE_DATA_OFFSET;
      layout->buffer_data_type[1] = NANOARROW_TYPE_INT64;
      layout->element_size_bits[1] = 64;
      layout->buffer_type[2] = NANOARROW_BUFFER_TYPE_DATA;
      layout->buffer_data_type[2] = NANOARROW_TYPE_STRING;
      break;
    case NANOARROW_TYPE_LARGE_BINARY:
      layout->buffer_type[1] = NANOARROW_BUFFER_TYPE_DATA_OFFSET;
      layout->buffer_data_type[1] = NANOARROW_TYPE_INT64;
      layout->element_size_bits[1] = 64;
      layout->buffer_type[2] = NANOARROW_BUFFER_TYPE_DATA;
      layout->buffer_data_type[2] = NANOARROW_TYPE_BINARY;
      break;

    default:
      break;
  }
}

void* ArrowMalloc(int64_t size) { return malloc(size); }

void* ArrowRealloc(void* ptr, int64_t size) { return realloc(ptr, size); }

void ArrowFree(void* ptr) { free(ptr); }

static uint8_t* ArrowBufferAllocatorMallocReallocate(
    struct ArrowBufferAllocator* allocator, uint8_t* ptr, int64_t old_size,
    int64_t new_size) {
  return (uint8_t*)ArrowRealloc(ptr, new_size);
}

static void ArrowBufferAllocatorMallocFree(struct ArrowBufferAllocator* allocator,
                                           uint8_t* ptr, int64_t size) {
  ArrowFree(ptr);
}

static struct ArrowBufferAllocator ArrowBufferAllocatorMalloc = {
    &ArrowBufferAllocatorMallocReallocate, &ArrowBufferAllocatorMallocFree, NULL};

struct ArrowBufferAllocator ArrowBufferAllocatorDefault(void) {
  return ArrowBufferAllocatorMalloc;
}

static uint8_t* ArrowBufferAllocatorNeverReallocate(
    struct ArrowBufferAllocator* allocator, uint8_t* ptr, int64_t old_size,
    int64_t new_size) {
  return NULL;
}

struct ArrowBufferAllocator ArrowBufferDeallocator(
    void (*custom_free)(struct ArrowBufferAllocator* allocator, uint8_t* ptr,
                        int64_t size),
    void* private_data) {
  struct ArrowBufferAllocator allocator;
  allocator.reallocate = &ArrowBufferAllocatorNeverReallocate;
  allocator.free = custom_free;
  allocator.private_data = private_data;
  return allocator;
}
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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nanoarrow.h"

static void ArrowSchemaRelease(struct ArrowSchema* schema) {
  if (schema->format != NULL) ArrowFree((void*)schema->format);
  if (schema->name != NULL) ArrowFree((void*)schema->name);
  if (schema->metadata != NULL) ArrowFree((void*)schema->metadata);

  // This object owns the memory for all the children, but those
  // children may have been generated elsewhere and might have
  // their own release() callback.
  if (schema->children != NULL) {
    for (int64_t i = 0; i < schema->n_children; i++) {
      if (schema->children[i] != NULL) {
        if (schema->children[i]->release != NULL) {
          schema->children[i]->release(schema->children[i]);
        }

        ArrowFree(schema->children[i]);
      }
    }

    ArrowFree(schema->children);
  }

  // This object owns the memory for the dictionary but it
  // may have been generated somewhere else and have its own
  // release() callback.
  if (schema->dictionary != NULL) {
    if (schema->dictionary->release != NULL) {
      schema->dictionary->release(schema->dictionary);
    }

    ArrowFree(schema->dictionary);
  }

  // private data not currently used
  if (schema->private_data != NULL) {
    ArrowFree(schema->private_data);
  }

  schema->release = NULL;
}

static const char* ArrowSchemaFormatTemplate(enum ArrowType type) {
  switch (type) {
    case NANOARROW_TYPE_UNINITIALIZED:
      return NULL;
    case NANOARROW_TYPE_NA:
      return "n";
    case NANOARROW_TYPE_BOOL:
      return "b";

    case NANOARROW_TYPE_UINT8:
      return "C";
    case NANOARROW_TYPE_INT8:
      return "c";
    case NANOARROW_TYPE_UINT16:
      return "S";
    case NANOARROW_TYPE_INT16:
      return "s";
    case NANOARROW_TYPE_UINT32:
      return "I";
    case NANOARROW_TYPE_INT32:
      return "i";
    case NANOARROW_TYPE_UINT64:
      return "L";
    case NANOARROW_TYPE_INT64:
      return "l";

    case NANOARROW_TYPE_HALF_FLOAT:
      return "e";
    case NANOARROW_TYPE_FLOAT:
      return "f";
    case NANOARROW_TYPE_DOUBLE:
      return "g";

    case NANOARROW_TYPE_STRING:
      return "u";
    case NANOARROW_TYPE_LARGE_STRING:
      return "U";
    case NANOARROW_TYPE_BINARY:
      return "z";
    case NANOARROW_TYPE_LARGE_BINARY:
      return "Z";

    case NANOARROW_TYPE_DATE32:
      return "tdD";
    case NANOARROW_TYPE_DATE64:
      return "tdm";
    case NANOARROW_TYPE_INTERVAL_MONTHS:
      return "tiM";
    case NANOARROW_TYPE_INTERVAL_DAY_TIME:
      return "tiD";
    case NANOARROW_TYPE_INTERVAL_MONTH_DAY_NANO:
      return "tin";

    case NANOARROW_TYPE_LIST:
      return "+l";
    case NANOARROW_TYPE_LARGE_LIST:
      return "+L";
    case NANOARROW_TYPE_STRUCT:
      return "+s";
    case NANOARROW_TYPE_MAP:
      return "+m";

    default:
      return NULL;
  }
}

static int ArrowSchemaInitChildrenIfNeeded(struct ArrowSchema* schema,
                                           enum ArrowType type) {
  switch (type) {
    case NANOARROW_TYPE_LIST:
    case NANOARROW_TYPE_LARGE_LIST:
    case NANOARROW_TYPE_FIXED_SIZE_LIST:
      NANOARROW_RETURN_NOT_OK(ArrowSchemaAllocateChildren(schema, 1));
      ArrowSchemaInit(schema->children[0]);
      NANOARROW_RETURN_NOT_OK(ArrowSchemaSetName(schema->children[0], "item"));
      break;
    case NANOARROW_TYPE_MAP:
      NANOARROW_RETURN_NOT_OK(ArrowSchemaAllocateChildren(schema, 1));
      NANOARROW_RETURN_NOT_OK(
          ArrowSchemaInitFromType(schema->children[0], NANOARROW_TYPE_STRUCT));
      NANOARROW_RETURN_NOT_OK(ArrowSchemaSetName(schema->children[0], "entries"));
      schema->children[0]->flags &= ~ARROW_FLAG_NULLABLE;
      NANOARROW_RETURN_NOT_OK(ArrowSchemaAllocateChildren(schema->children[0], 2));
      ArrowSchemaInit(schema->children[0]->children[0]);
      ArrowSchemaInit(schema->children[0]->children[1]);
      NANOARROW_RETURN_NOT_OK(
          ArrowSchemaSetName(schema->children[0]->children[0], "key"));
      schema->children[0]->children[0]->flags &= ~ARROW_FLAG_NULLABLE;
      NANOARROW_RETURN_NOT_OK(
          ArrowSchemaSetName(schema->children[0]->children[1], "value"));
      break;
    default:
      break;
  }

  return NANOARROW_OK;
}

void ArrowSchemaInit(struct ArrowSchema* schema) {
  schema->format = NULL;
  schema->name = NULL;
  schema->metadata = NULL;
  schema->flags = ARROW_FLAG_NULLABLE;
  schema->n_children = 0;
  schema->children = NULL;
  schema->dictionary = NULL;
  schema->private_data = NULL;
  schema->release = &ArrowSchemaRelease;
}

ArrowErrorCode ArrowSchemaSetType(struct ArrowSchema* schema, enum ArrowType type) {
  // We don't allocate the dictionary because it has to be nullptr
  // for non-dictionary-encoded arrays.

  // Set the format to a valid format string for type
  const char* template_format = ArrowSchemaFormatTemplate(type);

  // If type isn't recognized and not explicitly unset
  if (template_format == NULL && type != NANOARROW_TYPE_UNINITIALIZED) {
    return EINVAL;
  }

  NANOARROW_RETURN_NOT_OK(ArrowSchemaSetFormat(schema, template_format));

  // For types with an umabiguous child structure, allocate children
  return ArrowSchemaInitChildrenIfNeeded(schema, type);
}

ArrowErrorCode ArrowSchemaSetTypeStruct(struct ArrowSchema* schema, int64_t n_children) {
  NANOARROW_RETURN_NOT_OK(ArrowSchemaSetType(schema, NANOARROW_TYPE_STRUCT));
  NANOARROW_RETURN_NOT_OK(ArrowSchemaAllocateChildren(schema, n_children));
  for (int64_t i = 0; i < n_children; i++) {
    ArrowSchemaInit(schema->children[i]);
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowSchemaInitFromType(struct ArrowSchema* schema, enum ArrowType type) {
  ArrowSchemaInit(schema);

  int result = ArrowSchemaSetType(schema, type);
  if (result != NANOARROW_OK) {
    schema->release(schema);
    return result;
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowSchemaSetTypeFixedSize(struct ArrowSchema* schema,
                                           enum ArrowType type, int32_t fixed_size) {
  if (fixed_size <= 0) {
    return EINVAL;
  }

  char buffer[64];
  int n_chars;
  switch (type) {
    case NANOARROW_TYPE_FIXED_SIZE_BINARY:
      n_chars = snprintf(buffer, sizeof(buffer), "w:%d", (int)fixed_size);
      break;
    case NANOARROW_TYPE_FIXED_SIZE_LIST:
      n_chars = snprintf(buffer, sizeof(buffer), "+w:%d", (int)fixed_size);
      break;
    default:
      return EINVAL;
  }

  buffer[n_chars] = '\0';
  NANOARROW_RETURN_NOT_OK(ArrowSchemaSetFormat(schema, buffer));

  if (type == NANOARROW_TYPE_FIXED_SIZE_LIST) {
    NANOARROW_RETURN_NOT_OK(ArrowSchemaInitChildrenIfNeeded(schema, type));
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowSchemaSetTypeDecimal(struct ArrowSchema* schema, enum ArrowType type,
                                         int32_t decimal_precision,
                                         int32_t decimal_scale) {
  if (decimal_precision <= 0) {
    return EINVAL;
  }

  char buffer[64];
  int n_chars;
  switch (type) {
    case NANOARROW_TYPE_DECIMAL128:
      n_chars =
          snprintf(buffer, sizeof(buffer), "d:%d,%d", decimal_precision, decimal_scale);
      break;
    case NANOARROW_TYPE_DECIMAL256:
      n_chars = snprintf(buffer, sizeof(buffer), "d:%d,%d,256", decimal_precision,
                         decimal_scale);
      break;
    default:
      return EINVAL;
  }

  buffer[n_chars] = '\0';
  return ArrowSchemaSetFormat(schema, buffer);
}

static const char* ArrowTimeUnitFormatString(enum ArrowTimeUnit time_unit) {
  switch (time_unit) {
    case NANOARROW_TIME_UNIT_SECOND:
      return "s";
    case NANOARROW_TIME_UNIT_MILLI:
      return "m";
    case NANOARROW_TIME_UNIT_MICRO:
      return "u";
    case NANOARROW_TIME_UNIT_NANO:
      return "n";
    default:
      return NULL;
  }
}

ArrowErrorCode ArrowSchemaSetTypeDateTime(struct ArrowSchema* schema, enum ArrowType type,
                                          enum ArrowTimeUnit time_unit,
                                          const char* timezone) {
  const char* time_unit_str = ArrowTimeUnitFormatString(time_unit);
  if (time_unit_str == NULL) {
    return EINVAL;
  }

  char buffer[128];
  int n_chars;
  switch (type) {
    case NANOARROW_TYPE_TIME32:
    case NANOARROW_TYPE_TIME64:
      if (timezone != NULL) {
        return EINVAL;
      }
      n_chars = snprintf(buffer, sizeof(buffer), "tt%s", time_unit_str);
      break;
    case NANOARROW_TYPE_TIMESTAMP:
      if (timezone == NULL) {
        timezone = "";
      }
      n_chars = snprintf(buffer, sizeof(buffer), "ts%s:%s", time_unit_str, timezone);
      break;
    case NANOARROW_TYPE_DURATION:
      if (timezone != NULL) {
        return EINVAL;
      }
      n_chars = snprintf(buffer, sizeof(buffer), "tD%s", time_unit_str);
      break;
    default:
      return EINVAL;
  }

  if (((size_t)n_chars) >= sizeof(buffer)) {
    return ERANGE;
  }

  buffer[n_chars] = '\0';

  return ArrowSchemaSetFormat(schema, buffer);
}

ArrowErrorCode ArrowSchemaSetTypeUnion(struct ArrowSchema* schema, enum ArrowType type,
                                       int64_t n_children) {
  if (n_children < 0 || n_children > 127) {
    return EINVAL;
  }

  // Max valid size would be +ud:0,1,...126 = 401 characters + null terminator
  char format_out[512];
  int64_t format_out_size = 512;
  memset(format_out, 0, format_out_size);
  int n_chars;
  char* format_cursor = format_out;

  switch (type) {
    case NANOARROW_TYPE_SPARSE_UNION:
      n_chars = snprintf(format_cursor, format_out_size, "+us:");
      format_cursor += n_chars;
      format_out_size -= n_chars;
      break;
    case NANOARROW_TYPE_DENSE_UNION:
      n_chars = snprintf(format_cursor, format_out_size, "+ud:");
      format_cursor += n_chars;
      format_out_size -= n_chars;
      break;
    default:
      return EINVAL;
  }

  if (n_children > 0) {
    n_chars = snprintf(format_cursor, format_out_size, "0");
    format_cursor += n_chars;
    format_out_size -= n_chars;

    for (int64_t i = 1; i < n_children; i++) {
      n_chars = snprintf(format_cursor, format_out_size, ",%d", (int)i);
      format_cursor += n_chars;
      format_out_size -= n_chars;
    }
  }

  NANOARROW_RETURN_NOT_OK(ArrowSchemaSetFormat(schema, format_out));

  NANOARROW_RETURN_NOT_OK(ArrowSchemaAllocateChildren(schema, n_children));
  for (int64_t i = 0; i < n_children; i++) {
    ArrowSchemaInit(schema->children[i]);
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowSchemaSetFormat(struct ArrowSchema* schema, const char* format) {
  if (schema->format != NULL) {
    ArrowFree((void*)schema->format);
  }

  if (format != NULL) {
    size_t format_size = strlen(format) + 1;
    schema->format = (const char*)ArrowMalloc(format_size);
    if (schema->format == NULL) {
      return ENOMEM;
    }

    memcpy((void*)schema->format, format, format_size);
  } else {
    schema->format = NULL;
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowSchemaSetName(struct ArrowSchema* schema, const char* name) {
  if (schema->name != NULL) {
    ArrowFree((void*)schema->name);
  }

  if (name != NULL) {
    size_t name_size = strlen(name) + 1;
    schema->name = (const char*)ArrowMalloc(name_size);
    if (schema->name == NULL) {
      return ENOMEM;
    }

    memcpy((void*)schema->name, name, name_size);
  } else {
    schema->name = NULL;
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowSchemaSetMetadata(struct ArrowSchema* schema, const char* metadata) {
  if (schema->metadata != NULL) {
    ArrowFree((void*)schema->metadata);
  }

  if (metadata != NULL) {
    size_t metadata_size = ArrowMetadataSizeOf(metadata);
    schema->metadata = (const char*)ArrowMalloc(metadata_size);
    if (schema->metadata == NULL) {
      return ENOMEM;
    }

    memcpy((void*)schema->metadata, metadata, metadata_size);
  } else {
    schema->metadata = NULL;
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowSchemaAllocateChildren(struct ArrowSchema* schema,
                                           int64_t n_children) {
  if (schema->children != NULL) {
    return EEXIST;
  }

  if (n_children > 0) {
    schema->children =
        (struct ArrowSchema**)ArrowMalloc(n_children * sizeof(struct ArrowSchema*));

    if (schema->children == NULL) {
      return ENOMEM;
    }

    schema->n_children = n_children;

    memset(schema->children, 0, n_children * sizeof(struct ArrowSchema*));

    for (int64_t i = 0; i < n_children; i++) {
      schema->children[i] = (struct ArrowSchema*)ArrowMalloc(sizeof(struct ArrowSchema));

      if (schema->children[i] == NULL) {
        return ENOMEM;
      }

      schema->children[i]->release = NULL;
    }
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowSchemaAllocateDictionary(struct ArrowSchema* schema) {
  if (schema->dictionary != NULL) {
    return EEXIST;
  }

  schema->dictionary = (struct ArrowSchema*)ArrowMalloc(sizeof(struct ArrowSchema));
  if (schema->dictionary == NULL) {
    return ENOMEM;
  }

  schema->dictionary->release = NULL;
  return NANOARROW_OK;
}

ArrowErrorCode ArrowSchemaDeepCopy(struct ArrowSchema* schema,
                                   struct ArrowSchema* schema_out) {
  ArrowSchemaInit(schema_out);

  int result = ArrowSchemaSetFormat(schema_out, schema->format);
  if (result != NANOARROW_OK) {
    schema_out->release(schema_out);
    return result;
  }

  schema_out->flags = schema->flags;

  result = ArrowSchemaSetName(schema_out, schema->name);
  if (result != NANOARROW_OK) {
    schema_out->release(schema_out);
    return result;
  }

  result = ArrowSchemaSetMetadata(schema_out, schema->metadata);
  if (result != NANOARROW_OK) {
    schema_out->release(schema_out);
    return result;
  }

  result = ArrowSchemaAllocateChildren(schema_out, schema->n_children);
  if (result != NANOARROW_OK) {
    schema_out->release(schema_out);
    return result;
  }

  for (int64_t i = 0; i < schema->n_children; i++) {
    result = ArrowSchemaDeepCopy(schema->children[i], schema_out->children[i]);
    if (result != NANOARROW_OK) {
      schema_out->release(schema_out);
      return result;
    }
  }

  if (schema->dictionary != NULL) {
    result = ArrowSchemaAllocateDictionary(schema_out);
    if (result != NANOARROW_OK) {
      schema_out->release(schema_out);
      return result;
    }

    result = ArrowSchemaDeepCopy(schema->dictionary, schema_out->dictionary);
    if (result != NANOARROW_OK) {
      schema_out->release(schema_out);
      return result;
    }
  }

  return NANOARROW_OK;
}

static void ArrowSchemaViewSetPrimitive(struct ArrowSchemaView* schema_view,
                                        enum ArrowType type) {
  schema_view->type = type;
  schema_view->storage_type = type;
}

static ArrowErrorCode ArrowSchemaViewParse(struct ArrowSchemaView* schema_view,
                                           const char* format,
                                           const char** format_end_out,
                                           struct ArrowError* error) {
  *format_end_out = format;

  // needed for decimal parsing
  const char* parse_start;
  char* parse_end;

  switch (format[0]) {
    case 'n':
      schema_view->type = NANOARROW_TYPE_NA;
      schema_view->storage_type = NANOARROW_TYPE_NA;
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 'b':
      ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_BOOL);
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 'c':
      ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT8);
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 'C':
      ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_UINT8);
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 's':
      ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT16);
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 'S':
      ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_UINT16);
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 'i':
      ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT32);
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 'I':
      ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_UINT32);
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 'l':
      ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT64);
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 'L':
      ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_UINT64);
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 'e':
      ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_HALF_FLOAT);
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 'f':
      ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_FLOAT);
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 'g':
      ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_DOUBLE);
      *format_end_out = format + 1;
      return NANOARROW_OK;

    // decimal
    case 'd':
      if (format[1] != ':' || format[2] == '\0') {
        ArrowErrorSet(error, "Expected ':precision,scale[,bitwidth]' following 'd'",
                      format + 3);
        return EINVAL;
      }

      parse_start = format + 2;
      schema_view->decimal_precision = (int32_t)strtol(parse_start, &parse_end, 10);
      if (parse_end == parse_start || parse_end[0] != ',') {
        ArrowErrorSet(error, "Expected 'precision,scale[,bitwidth]' following 'd:'");
        return EINVAL;
      }

      parse_start = parse_end + 1;
      schema_view->decimal_scale = (int32_t)strtol(parse_start, &parse_end, 10);
      if (parse_end == parse_start) {
        ArrowErrorSet(error, "Expected 'scale[,bitwidth]' following 'd:precision,'");
        return EINVAL;
      } else if (parse_end[0] != ',') {
        schema_view->decimal_bitwidth = 128;
      } else {
        parse_start = parse_end + 1;
        schema_view->decimal_bitwidth = (int32_t)strtol(parse_start, &parse_end, 10);
        if (parse_start == parse_end) {
          ArrowErrorSet(error, "Expected precision following 'd:precision,scale,'");
          return EINVAL;
        }
      }

      *format_end_out = parse_end;

      switch (schema_view->decimal_bitwidth) {
        case 128:
          ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_DECIMAL128);
          return NANOARROW_OK;
        case 256:
          ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_DECIMAL256);
          return NANOARROW_OK;
        default:
          ArrowErrorSet(error, "Expected decimal bitwidth of 128 or 256 but found %d",
                        (int)schema_view->decimal_bitwidth);
          return EINVAL;
      }

    // validity + data
    case 'w':
      schema_view->type = NANOARROW_TYPE_FIXED_SIZE_BINARY;
      schema_view->storage_type = NANOARROW_TYPE_FIXED_SIZE_BINARY;
      if (format[1] != ':' || format[2] == '\0') {
        ArrowErrorSet(error, "Expected ':<width>' following 'w'");
        return EINVAL;
      }

      schema_view->fixed_size = (int32_t)strtol(format + 2, (char**)format_end_out, 10);
      return NANOARROW_OK;

    // validity + offset + data
    case 'z':
      schema_view->type = NANOARROW_TYPE_BINARY;
      schema_view->storage_type = NANOARROW_TYPE_BINARY;
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 'u':
      schema_view->type = NANOARROW_TYPE_STRING;
      schema_view->storage_type = NANOARROW_TYPE_STRING;
      *format_end_out = format + 1;
      return NANOARROW_OK;

    // validity + large_offset + data
    case 'Z':
      schema_view->type = NANOARROW_TYPE_LARGE_BINARY;
      schema_view->storage_type = NANOARROW_TYPE_LARGE_BINARY;
      *format_end_out = format + 1;
      return NANOARROW_OK;
    case 'U':
      schema_view->type = NANOARROW_TYPE_LARGE_STRING;
      schema_view->storage_type = NANOARROW_TYPE_LARGE_STRING;
      *format_end_out = format + 1;
      return NANOARROW_OK;

    // nested types
    case '+':
      switch (format[1]) {
        // list has validity + offset or offset
        case 'l':
          schema_view->storage_type = NANOARROW_TYPE_LIST;
          schema_view->type = NANOARROW_TYPE_LIST;
          *format_end_out = format + 2;
          return NANOARROW_OK;

        // large list has validity + large_offset or large_offset
        case 'L':
          schema_view->storage_type = NANOARROW_TYPE_LARGE_LIST;
          schema_view->type = NANOARROW_TYPE_LARGE_LIST;
          *format_end_out = format + 2;
          return NANOARROW_OK;

        // just validity buffer
        case 'w':
          if (format[2] != ':' || format[3] == '\0') {
            ArrowErrorSet(error, "Expected ':<width>' following '+w'");
            return EINVAL;
          }

          schema_view->storage_type = NANOARROW_TYPE_FIXED_SIZE_LIST;
          schema_view->type = NANOARROW_TYPE_FIXED_SIZE_LIST;
          schema_view->fixed_size =
              (int32_t)strtol(format + 3, (char**)format_end_out, 10);
          return NANOARROW_OK;
        case 's':
          schema_view->storage_type = NANOARROW_TYPE_STRUCT;
          schema_view->type = NANOARROW_TYPE_STRUCT;
          *format_end_out = format + 2;
          return NANOARROW_OK;
        case 'm':
          schema_view->storage_type = NANOARROW_TYPE_MAP;
          schema_view->type = NANOARROW_TYPE_MAP;
          *format_end_out = format + 2;
          return NANOARROW_OK;

        // unions
        case 'u':
          switch (format[2]) {
            case 'd':
              schema_view->storage_type = NANOARROW_TYPE_DENSE_UNION;
              schema_view->type = NANOARROW_TYPE_DENSE_UNION;
              break;
            case 's':
              schema_view->storage_type = NANOARROW_TYPE_SPARSE_UNION;
              schema_view->type = NANOARROW_TYPE_SPARSE_UNION;
              break;
            default:
              ArrowErrorSet(error,
                            "Expected union format string +us:<type_ids> or "
                            "+ud:<type_ids> but found '%s'",
                            format);
              return EINVAL;
          }

          if (format[3] == ':') {
            schema_view->union_type_ids = format + 4;
            int64_t n_type_ids =
                _ArrowParseUnionTypeIds(schema_view->union_type_ids, NULL);
            if (n_type_ids != schema_view->schema->n_children) {
              ArrowErrorSet(
                  error,
                  "Expected union type_ids parameter to be a comma-separated list of %ld "
                  "values between 0 and 127 but found '%s'",
                  (long)schema_view->schema->n_children, schema_view->union_type_ids);
              return EINVAL;
            }
            *format_end_out = format + strlen(format);
            return NANOARROW_OK;
          } else {
            ArrowErrorSet(error,
                          "Expected union format string +us:<type_ids> or +ud:<type_ids> "
                          "but found '%s'",
                          format);
            return EINVAL;
          }

        default:
          ArrowErrorSet(error, "Expected nested type format string but found '%s'",
                        format);
          return EINVAL;
      }

    // date/time types
    case 't':
      switch (format[1]) {
        // date
        case 'd':
          switch (format[2]) {
            case 'D':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT32);
              schema_view->type = NANOARROW_TYPE_DATE32;
              *format_end_out = format + 3;
              return NANOARROW_OK;
            case 'm':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT64);
              schema_view->type = NANOARROW_TYPE_DATE64;
              *format_end_out = format + 3;
              return NANOARROW_OK;
            default:
              ArrowErrorSet(error, "Expected 'D' or 'm' following 'td' but found '%s'",
                            format + 2);
              return EINVAL;
          }

        // time of day
        case 't':
          switch (format[2]) {
            case 's':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT32);
              schema_view->type = NANOARROW_TYPE_TIME32;
              schema_view->time_unit = NANOARROW_TIME_UNIT_SECOND;
              *format_end_out = format + 3;
              return NANOARROW_OK;
            case 'm':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT32);
              schema_view->type = NANOARROW_TYPE_TIME32;
              schema_view->time_unit = NANOARROW_TIME_UNIT_MILLI;
              *format_end_out = format + 3;
              return NANOARROW_OK;
            case 'u':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT64);
              schema_view->type = NANOARROW_TYPE_TIME64;
              schema_view->time_unit = NANOARROW_TIME_UNIT_MICRO;
              *format_end_out = format + 3;
              return NANOARROW_OK;
            case 'n':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT64);
              schema_view->type = NANOARROW_TYPE_TIME64;
              schema_view->time_unit = NANOARROW_TIME_UNIT_NANO;
              *format_end_out = format + 3;
              return NANOARROW_OK;
            default:
              ArrowErrorSet(
                  error, "Expected 's', 'm', 'u', or 'n' following 'tt' but found '%s'",
                  format + 2);
              return EINVAL;
          }

        // timestamp
        case 's':
          switch (format[2]) {
            case 's':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT64);
              schema_view->type = NANOARROW_TYPE_TIMESTAMP;
              schema_view->time_unit = NANOARROW_TIME_UNIT_SECOND;
              break;
            case 'm':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT64);
              schema_view->type = NANOARROW_TYPE_TIMESTAMP;
              schema_view->time_unit = NANOARROW_TIME_UNIT_MILLI;
              break;
            case 'u':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT64);
              schema_view->type = NANOARROW_TYPE_TIMESTAMP;
              schema_view->time_unit = NANOARROW_TIME_UNIT_MICRO;
              break;
            case 'n':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT64);
              schema_view->type = NANOARROW_TYPE_TIMESTAMP;
              schema_view->time_unit = NANOARROW_TIME_UNIT_NANO;
              break;
            default:
              ArrowErrorSet(
                  error, "Expected 's', 'm', 'u', or 'n' following 'ts' but found '%s'",
                  format + 2);
              return EINVAL;
          }

          if (format[3] != ':') {
            ArrowErrorSet(error, "Expected ':' following '%.3s' but found '%s'", format,
                          format + 3);
            return EINVAL;
          }

          schema_view->timezone = format + 4;
          *format_end_out = format + strlen(format);
          return NANOARROW_OK;

        // duration
        case 'D':
          switch (format[2]) {
            case 's':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT64);
              schema_view->type = NANOARROW_TYPE_DURATION;
              schema_view->time_unit = NANOARROW_TIME_UNIT_SECOND;
              *format_end_out = format + 3;
              return NANOARROW_OK;
            case 'm':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT64);
              schema_view->type = NANOARROW_TYPE_DURATION;
              schema_view->time_unit = NANOARROW_TIME_UNIT_MILLI;
              *format_end_out = format + 3;
              return NANOARROW_OK;
            case 'u':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT64);
              schema_view->type = NANOARROW_TYPE_DURATION;
              schema_view->time_unit = NANOARROW_TIME_UNIT_MICRO;
              *format_end_out = format + 3;
              return NANOARROW_OK;
            case 'n':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INT64);
              schema_view->type = NANOARROW_TYPE_DURATION;
              schema_view->time_unit = NANOARROW_TIME_UNIT_NANO;
              *format_end_out = format + 3;
              return NANOARROW_OK;
            default:
              ArrowErrorSet(error,
                            "Expected 's', 'm', u', or 'n' following 'tD' but found '%s'",
                            format + 2);
              return EINVAL;
          }

        // interval
        case 'i':
          switch (format[2]) {
            case 'M':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INTERVAL_MONTHS);
              *format_end_out = format + 3;
              return NANOARROW_OK;
            case 'D':
              ArrowSchemaViewSetPrimitive(schema_view, NANOARROW_TYPE_INTERVAL_DAY_TIME);
              *format_end_out = format + 3;
              return NANOARROW_OK;
            case 'n':
              ArrowSchemaViewSetPrimitive(schema_view,
                                          NANOARROW_TYPE_INTERVAL_MONTH_DAY_NANO);
              *format_end_out = format + 3;
              return NANOARROW_OK;
            default:
              ArrowErrorSet(error,
                            "Expected 'M', 'D', or 'n' following 'ti' but found '%s'",
                            format + 2);
              return EINVAL;
          }

        default:
          ArrowErrorSet(
              error, "Expected 'd', 't', 's', 'D', or 'i' following 't' but found '%s'",
              format + 1);
          return EINVAL;
      }

    default:
      ArrowErrorSet(error, "Unknown format: '%s'", format);
      return EINVAL;
  }
}

static ArrowErrorCode ArrowSchemaViewValidateNChildren(
    struct ArrowSchemaView* schema_view, int64_t n_children, struct ArrowError* error) {
  if (n_children != -1 && schema_view->schema->n_children != n_children) {
    ArrowErrorSet(error, "Expected schema with %d children but found %d children",
                  (int)n_children, (int)schema_view->schema->n_children);
    return EINVAL;
  }

  // Don't do a full validation of children but do check that they won't
  // segfault if inspected
  struct ArrowSchema* child;
  for (int64_t i = 0; i < schema_view->schema->n_children; i++) {
    child = schema_view->schema->children[i];
    if (child == NULL) {
      ArrowErrorSet(error, "Expected valid schema at schema->children[%d] but found NULL",
                    i);
      return EINVAL;
    } else if (child->release == NULL) {
      ArrowErrorSet(
          error,
          "Expected valid schema at schema->children[%d] but found a released schema", i);
      return EINVAL;
    }
  }

  return NANOARROW_OK;
}

static ArrowErrorCode ArrowSchemaViewValidateUnion(struct ArrowSchemaView* schema_view,
                                                   struct ArrowError* error) {
  return ArrowSchemaViewValidateNChildren(schema_view, -1, error);
}

static ArrowErrorCode ArrowSchemaViewValidateMap(struct ArrowSchemaView* schema_view,
                                                 struct ArrowError* error) {
  NANOARROW_RETURN_NOT_OK(ArrowSchemaViewValidateNChildren(schema_view, 1, error));

  if (schema_view->schema->children[0]->n_children != 2) {
    ArrowErrorSet(error, "Expected child of map type to have 2 children but found %d",
                  (int)schema_view->schema->children[0]->n_children);
    return EINVAL;
  }

  if (strcmp(schema_view->schema->children[0]->format, "+s") != 0) {
    ArrowErrorSet(error, "Expected format of child of map type to be '+s' but found '%s'",
                  schema_view->schema->children[0]->format);
    return EINVAL;
  }

  if (schema_view->schema->children[0]->flags & ARROW_FLAG_NULLABLE) {
    ArrowErrorSet(error,
                  "Expected child of map type to be non-nullable but was nullable");
    return EINVAL;
  }

  if (schema_view->schema->children[0]->children[0]->flags & ARROW_FLAG_NULLABLE) {
    ArrowErrorSet(error, "Expected key of map type to be non-nullable but was nullable");
    return EINVAL;
  }

  return NANOARROW_OK;
}

static ArrowErrorCode ArrowSchemaViewValidateDictionary(
    struct ArrowSchemaView* schema_view, struct ArrowError* error) {
  // check for valid index type
  switch (schema_view->storage_type) {
    case NANOARROW_TYPE_UINT8:
    case NANOARROW_TYPE_INT8:
    case NANOARROW_TYPE_UINT16:
    case NANOARROW_TYPE_INT16:
    case NANOARROW_TYPE_UINT32:
    case NANOARROW_TYPE_INT32:
    case NANOARROW_TYPE_UINT64:
    case NANOARROW_TYPE_INT64:
      break;
    default:
      ArrowErrorSet(
          error,
          "Expected dictionary schema index type to be an integral type but found '%s'",
          schema_view->schema->format);
      return EINVAL;
  }

  struct ArrowSchemaView dictionary_schema_view;
  return ArrowSchemaViewInit(&dictionary_schema_view, schema_view->schema->dictionary,
                             error);
}

static ArrowErrorCode ArrowSchemaViewValidate(struct ArrowSchemaView* schema_view,
                                              enum ArrowType type,
                                              struct ArrowError* error) {
  switch (type) {
    case NANOARROW_TYPE_NA:
    case NANOARROW_TYPE_BOOL:
    case NANOARROW_TYPE_UINT8:
    case NANOARROW_TYPE_INT8:
    case NANOARROW_TYPE_UINT16:
    case NANOARROW_TYPE_INT16:
    case NANOARROW_TYPE_UINT32:
    case NANOARROW_TYPE_INT32:
    case NANOARROW_TYPE_UINT64:
    case NANOARROW_TYPE_INT64:
    case NANOARROW_TYPE_HALF_FLOAT:
    case NANOARROW_TYPE_FLOAT:
    case NANOARROW_TYPE_DOUBLE:
    case NANOARROW_TYPE_DECIMAL128:
    case NANOARROW_TYPE_DECIMAL256:
    case NANOARROW_TYPE_STRING:
    case NANOARROW_TYPE_LARGE_STRING:
    case NANOARROW_TYPE_BINARY:
    case NANOARROW_TYPE_LARGE_BINARY:
    case NANOARROW_TYPE_DATE32:
    case NANOARROW_TYPE_DATE64:
    case NANOARROW_TYPE_INTERVAL_MONTHS:
    case NANOARROW_TYPE_INTERVAL_DAY_TIME:
    case NANOARROW_TYPE_INTERVAL_MONTH_DAY_NANO:
    case NANOARROW_TYPE_TIMESTAMP:
    case NANOARROW_TYPE_TIME32:
    case NANOARROW_TYPE_TIME64:
    case NANOARROW_TYPE_DURATION:
      return ArrowSchemaViewValidateNChildren(schema_view, 0, error);

    case NANOARROW_TYPE_FIXED_SIZE_BINARY:
      if (schema_view->fixed_size <= 0) {
        ArrowErrorSet(error, "Expected size > 0 for fixed size binary but found size %d",
                      schema_view->fixed_size);
        return EINVAL;
      }
      return ArrowSchemaViewValidateNChildren(schema_view, 0, error);

    case NANOARROW_TYPE_LIST:
    case NANOARROW_TYPE_LARGE_LIST:
    case NANOARROW_TYPE_FIXED_SIZE_LIST:
      return ArrowSchemaViewValidateNChildren(schema_view, 1, error);

    case NANOARROW_TYPE_STRUCT:
      return ArrowSchemaViewValidateNChildren(schema_view, -1, error);

    case NANOARROW_TYPE_SPARSE_UNION:
    case NANOARROW_TYPE_DENSE_UNION:
      return ArrowSchemaViewValidateUnion(schema_view, error);

    case NANOARROW_TYPE_MAP:
      return ArrowSchemaViewValidateMap(schema_view, error);

    case NANOARROW_TYPE_DICTIONARY:
      return ArrowSchemaViewValidateDictionary(schema_view, error);

    default:
      ArrowErrorSet(error, "Expected a valid enum ArrowType value but found %d",
                    (int)schema_view->type);
      return EINVAL;
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowSchemaViewInit(struct ArrowSchemaView* schema_view,
                                   struct ArrowSchema* schema, struct ArrowError* error) {
  if (schema == NULL) {
    ArrowErrorSet(error, "Expected non-NULL schema");
    return EINVAL;
  }

  if (schema->release == NULL) {
    ArrowErrorSet(error, "Expected non-released schema");
    return EINVAL;
  }

  schema_view->schema = schema;

  const char* format = schema->format;
  if (format == NULL) {
    ArrowErrorSet(
        error,
        "Error parsing schema->format: Expected a null-terminated string but found NULL");
    return EINVAL;
  }

  size_t format_len = strlen(format);
  if (format_len == 0) {
    ArrowErrorSet(error, "Error parsing schema->format: Expected a string with size > 0");
    return EINVAL;
  }

  const char* format_end_out;
  ArrowErrorCode result =
      ArrowSchemaViewParse(schema_view, format, &format_end_out, error);

  if (result != NANOARROW_OK) {
    if (error != NULL) {
      char child_error[1024];
      memcpy(child_error, ArrowErrorMessage(error), 1024);
      ArrowErrorSet(error, "Error parsing schema->format: %s", child_error);
    }

    return result;
  }

  if ((format + format_len) != format_end_out) {
    ArrowErrorSet(error, "Error parsing schema->format '%s': parsed %d/%d characters",
                  format, (int)(format_end_out - format), (int)(format_len));
    return EINVAL;
  }

  if (schema->dictionary != NULL) {
    schema_view->type = NANOARROW_TYPE_DICTIONARY;
  }

  result = ArrowSchemaViewValidate(schema_view, schema_view->storage_type, error);
  if (result != NANOARROW_OK) {
    return result;
  }

  if (schema_view->storage_type != schema_view->type) {
    result = ArrowSchemaViewValidate(schema_view, schema_view->type, error);
    if (result != NANOARROW_OK) {
      return result;
    }
  }

  ArrowLayoutInit(&schema_view->layout, schema_view->storage_type);
  if (schema_view->storage_type == NANOARROW_TYPE_FIXED_SIZE_BINARY) {
    schema_view->layout.element_size_bits[1] = schema_view->fixed_size * 8;
  } else if (schema_view->storage_type == NANOARROW_TYPE_FIXED_SIZE_LIST) {
    schema_view->layout.child_size_elements = schema_view->fixed_size;
  }

  schema_view->extension_name = ArrowCharView(NULL);
  schema_view->extension_metadata = ArrowCharView(NULL);
  ArrowMetadataGetValue(schema->metadata, ArrowCharView("ARROW:extension:name"),
                        &schema_view->extension_name);
  ArrowMetadataGetValue(schema->metadata, ArrowCharView("ARROW:extension:metadata"),
                        &schema_view->extension_metadata);

  return NANOARROW_OK;
}

static int64_t ArrowSchemaTypeToStringInternal(struct ArrowSchemaView* schema_view,
                                               char* out, int64_t n) {
  const char* type_string = ArrowTypeString(schema_view->type);
  switch (schema_view->type) {
    case NANOARROW_TYPE_DECIMAL128:
    case NANOARROW_TYPE_DECIMAL256:
      return snprintf(out, n, "%s(%d, %d)", type_string,
                      (int)schema_view->decimal_precision,
                      (int)schema_view->decimal_scale);
    case NANOARROW_TYPE_TIMESTAMP:
      return snprintf(out, n, "%s('%s', '%s')", type_string,
                      ArrowTimeUnitString(schema_view->time_unit), schema_view->timezone);
    case NANOARROW_TYPE_TIME32:
    case NANOARROW_TYPE_TIME64:
    case NANOARROW_TYPE_DURATION:
      return snprintf(out, n, "%s('%s')", type_string,
                      ArrowTimeUnitString(schema_view->time_unit));
    case NANOARROW_TYPE_FIXED_SIZE_BINARY:
    case NANOARROW_TYPE_FIXED_SIZE_LIST:
      return snprintf(out, n, "%s(%ld)", type_string, (long)schema_view->fixed_size);
    case NANOARROW_TYPE_SPARSE_UNION:
    case NANOARROW_TYPE_DENSE_UNION:
      return snprintf(out, n, "%s([%s])", type_string, schema_view->union_type_ids);
    default:
      return snprintf(out, n, "%s", type_string);
  }
}

// Helper for bookkeeping to emulate sprintf()-like behaviour spread
// among multiple sprintf calls.
static inline void ArrowToStringLogChars(char** out, int64_t n_chars_last,
                                         int64_t* n_remaining, int64_t* n_chars) {
  *n_chars += n_chars_last;
  *n_remaining -= n_chars_last;

  // n_remaining is never less than 0
  if (*n_remaining < 0) {
    *n_remaining = 0;
  }

  // Can't do math on a NULL pointer
  if (*out != NULL) {
    *out += n_chars_last;
  }
}

int64_t ArrowSchemaToString(struct ArrowSchema* schema, char* out, int64_t n,
                            char recursive) {
  if (schema == NULL) {
    return snprintf(out, n, "[invalid: pointer is null]");
  }

  if (schema->release == NULL) {
    return snprintf(out, n, "[invalid: schema is released]");
  }

  struct ArrowSchemaView schema_view;
  struct ArrowError error;

  if (ArrowSchemaViewInit(&schema_view, schema, &error) != NANOARROW_OK) {
    return snprintf(out, n, "[invalid: %s]", ArrowErrorMessage(&error));
  }

  // Extension type and dictionary should include both the top-level type
  // and the storage type.
  int is_extension = schema_view.extension_name.size_bytes > 0;
  int is_dictionary = schema->dictionary != NULL;
  int64_t n_chars = 0;
  int64_t n_chars_last = 0;

  // Uncommon but not technically impossible that both are true
  if (is_extension && is_dictionary) {
    n_chars_last = snprintf(
        out, n, "%.*s{dictionary(%s)<", (int)schema_view.extension_name.size_bytes,
        schema_view.extension_name.data, ArrowTypeString(schema_view.storage_type));
  } else if (is_extension) {
    n_chars_last = snprintf(out, n, "%.*s{", (int)schema_view.extension_name.size_bytes,
                            schema_view.extension_name.data);
  } else if (is_dictionary) {
    n_chars_last =
        snprintf(out, n, "dictionary(%s)<", ArrowTypeString(schema_view.storage_type));
  }

  ArrowToStringLogChars(&out, n_chars_last, &n, &n_chars);

  if (!is_dictionary) {
    n_chars_last = ArrowSchemaTypeToStringInternal(&schema_view, out, n);
  } else {
    n_chars_last = ArrowSchemaToString(schema->dictionary, out, n, recursive);
  }

  ArrowToStringLogChars(&out, n_chars_last, &n, &n_chars);

  if (recursive && schema->format[0] == '+') {
    n_chars_last = snprintf(out, n, "<");
    ArrowToStringLogChars(&out, n_chars_last, &n, &n_chars);

    for (int64_t i = 0; i < schema->n_children; i++) {
      if (i > 0) {
        n_chars_last = snprintf(out, n, ", ");
        ArrowToStringLogChars(&out, n_chars_last, &n, &n_chars);
      }

      // ArrowSchemaToStringInternal() will validate the child and print the error,
      // but we need the name first
      if (schema->children[i] != NULL && schema->children[i]->release != NULL &&
          schema->children[i]->name != NULL) {
        n_chars_last = snprintf(out, n, "%s: ", schema->children[i]->name);
        ArrowToStringLogChars(&out, n_chars_last, &n, &n_chars);
      }

      n_chars_last = ArrowSchemaToString(schema->children[i], out, n, recursive);
      ArrowToStringLogChars(&out, n_chars_last, &n, &n_chars);
    }

    n_chars_last = snprintf(out, n, ">");
    ArrowToStringLogChars(&out, n_chars_last, &n, &n_chars);
  }

  if (is_extension && is_dictionary) {
    n_chars += snprintf(out, n, ">}");
  } else if (is_extension) {
    n_chars += snprintf(out, n, "}");
  } else if (is_dictionary) {
    n_chars += snprintf(out, n, ">");
  }

  return n_chars;
}

ArrowErrorCode ArrowMetadataReaderInit(struct ArrowMetadataReader* reader,
                                       const char* metadata) {
  reader->metadata = metadata;

  if (reader->metadata == NULL) {
    reader->offset = 0;
    reader->remaining_keys = 0;
  } else {
    memcpy(&reader->remaining_keys, reader->metadata, sizeof(int32_t));
    reader->offset = sizeof(int32_t);
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowMetadataReaderRead(struct ArrowMetadataReader* reader,
                                       struct ArrowStringView* key_out,
                                       struct ArrowStringView* value_out) {
  if (reader->remaining_keys <= 0) {
    return EINVAL;
  }

  int64_t pos = 0;

  int32_t key_size;
  memcpy(&key_size, reader->metadata + reader->offset + pos, sizeof(int32_t));
  pos += sizeof(int32_t);

  key_out->data = reader->metadata + reader->offset + pos;
  key_out->size_bytes = key_size;
  pos += key_size;

  int32_t value_size;
  memcpy(&value_size, reader->metadata + reader->offset + pos, sizeof(int32_t));
  pos += sizeof(int32_t);

  value_out->data = reader->metadata + reader->offset + pos;
  value_out->size_bytes = value_size;
  pos += value_size;

  reader->offset += pos;
  reader->remaining_keys--;
  return NANOARROW_OK;
}

int64_t ArrowMetadataSizeOf(const char* metadata) {
  if (metadata == NULL) {
    return 0;
  }

  struct ArrowMetadataReader reader;
  struct ArrowStringView key;
  struct ArrowStringView value;
  ArrowMetadataReaderInit(&reader, metadata);

  int64_t size = sizeof(int32_t);
  while (ArrowMetadataReaderRead(&reader, &key, &value) == NANOARROW_OK) {
    size += sizeof(int32_t) + key.size_bytes + sizeof(int32_t) + value.size_bytes;
  }

  return size;
}

static ArrowErrorCode ArrowMetadataGetValueInternal(const char* metadata,
                                                    struct ArrowStringView* key,
                                                    struct ArrowStringView* value_out) {
  struct ArrowMetadataReader reader;
  struct ArrowStringView existing_key;
  struct ArrowStringView existing_value;
  ArrowMetadataReaderInit(&reader, metadata);

  while (ArrowMetadataReaderRead(&reader, &existing_key, &existing_value) ==
         NANOARROW_OK) {
    int key_equal = key->size_bytes == existing_key.size_bytes &&
                    strncmp(key->data, existing_key.data, existing_key.size_bytes) == 0;
    if (key_equal) {
      value_out->data = existing_value.data;
      value_out->size_bytes = existing_value.size_bytes;
      break;
    }
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowMetadataGetValue(const char* metadata, struct ArrowStringView key,
                                     struct ArrowStringView* value_out) {
  if (value_out == NULL) {
    return EINVAL;
  }

  return ArrowMetadataGetValueInternal(metadata, &key, value_out);
}

char ArrowMetadataHasKey(const char* metadata, struct ArrowStringView key) {
  struct ArrowStringView value = ArrowCharView(NULL);
  ArrowMetadataGetValue(metadata, key, &value);
  return value.data != NULL;
}

ArrowErrorCode ArrowMetadataBuilderInit(struct ArrowBuffer* buffer,
                                        const char* metadata) {
  ArrowBufferInit(buffer);
  return ArrowBufferAppend(buffer, metadata, ArrowMetadataSizeOf(metadata));
}

static ArrowErrorCode ArrowMetadataBuilderAppendInternal(struct ArrowBuffer* buffer,
                                                         struct ArrowStringView* key,
                                                         struct ArrowStringView* value) {
  if (value == NULL) {
    return NANOARROW_OK;
  }

  if (buffer->capacity_bytes == 0) {
    NANOARROW_RETURN_NOT_OK(ArrowBufferAppendInt32(buffer, 0));
  }

  if (((size_t)buffer->capacity_bytes) < sizeof(int32_t)) {
    return EINVAL;
  }

  int32_t n_keys;
  memcpy(&n_keys, buffer->data, sizeof(int32_t));

  int32_t key_size = (int32_t)key->size_bytes;
  int32_t value_size = (int32_t)value->size_bytes;
  NANOARROW_RETURN_NOT_OK(ArrowBufferReserve(
      buffer, sizeof(int32_t) + key_size + sizeof(int32_t) + value_size));

  ArrowBufferAppendUnsafe(buffer, &key_size, sizeof(int32_t));
  ArrowBufferAppendUnsafe(buffer, key->data, key_size);
  ArrowBufferAppendUnsafe(buffer, &value_size, sizeof(int32_t));
  ArrowBufferAppendUnsafe(buffer, value->data, value_size);

  n_keys++;
  memcpy(buffer->data, &n_keys, sizeof(int32_t));

  return NANOARROW_OK;
}

static ArrowErrorCode ArrowMetadataBuilderSetInternal(struct ArrowBuffer* buffer,
                                                      struct ArrowStringView* key,
                                                      struct ArrowStringView* value) {
  // Inspect the current value to see if we can avoid copying the buffer
  struct ArrowStringView current_value = ArrowCharView(NULL);
  NANOARROW_RETURN_NOT_OK(
      ArrowMetadataGetValueInternal((const char*)buffer->data, key, &current_value));

  // The key should be removed but no key exists
  if (value == NULL && current_value.data == NULL) {
    return NANOARROW_OK;
  }

  // The key/value can be appended because no key exists
  if (value != NULL && current_value.data == NULL) {
    return ArrowMetadataBuilderAppendInternal(buffer, key, value);
  }

  struct ArrowMetadataReader reader;
  struct ArrowStringView existing_key;
  struct ArrowStringView existing_value;
  NANOARROW_RETURN_NOT_OK(ArrowMetadataReaderInit(&reader, (const char*)buffer->data));

  struct ArrowBuffer new_buffer;
  NANOARROW_RETURN_NOT_OK(ArrowMetadataBuilderInit(&new_buffer, NULL));

  while (reader.remaining_keys > 0) {
    int result = ArrowMetadataReaderRead(&reader, &existing_key, &existing_value);
    if (result != NANOARROW_OK) {
      ArrowBufferReset(&new_buffer);
      return result;
    }

    if (key->size_bytes == existing_key.size_bytes &&
        strncmp((const char*)key->data, (const char*)existing_key.data,
                existing_key.size_bytes) == 0) {
      result = ArrowMetadataBuilderAppendInternal(&new_buffer, key, value);
      value = NULL;
    } else {
      result =
          ArrowMetadataBuilderAppendInternal(&new_buffer, &existing_key, &existing_value);
    }

    if (result != NANOARROW_OK) {
      ArrowBufferReset(&new_buffer);
      return result;
    }
  }

  ArrowBufferReset(buffer);
  ArrowBufferMove(&new_buffer, buffer);
  return NANOARROW_OK;
}

ArrowErrorCode ArrowMetadataBuilderAppend(struct ArrowBuffer* buffer,
                                          struct ArrowStringView key,
                                          struct ArrowStringView value) {
  return ArrowMetadataBuilderAppendInternal(buffer, &key, &value);
}

ArrowErrorCode ArrowMetadataBuilderSet(struct ArrowBuffer* buffer,
                                       struct ArrowStringView key,
                                       struct ArrowStringView value) {
  return ArrowMetadataBuilderSetInternal(buffer, &key, &value);
}

ArrowErrorCode ArrowMetadataBuilderRemove(struct ArrowBuffer* buffer,
                                          struct ArrowStringView key) {
  return ArrowMetadataBuilderSetInternal(buffer, &key, NULL);
}
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

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "nanoarrow.h"

static void ArrowArrayRelease(struct ArrowArray* array) {
  // Release buffers held by this array
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;
  if (private_data != NULL) {
    ArrowBitmapReset(&private_data->bitmap);
    ArrowBufferReset(&private_data->buffers[0]);
    ArrowBufferReset(&private_data->buffers[1]);
    ArrowFree(private_data);
  }

  // This object owns the memory for all the children, but those
  // children may have been generated elsewhere and might have
  // their own release() callback.
  if (array->children != NULL) {
    for (int64_t i = 0; i < array->n_children; i++) {
      if (array->children[i] != NULL) {
        if (array->children[i]->release != NULL) {
          array->children[i]->release(array->children[i]);
        }

        ArrowFree(array->children[i]);
      }
    }

    ArrowFree(array->children);
  }

  // This object owns the memory for the dictionary but it
  // may have been generated somewhere else and have its own
  // release() callback.
  if (array->dictionary != NULL) {
    if (array->dictionary->release != NULL) {
      array->dictionary->release(array->dictionary);
    }

    ArrowFree(array->dictionary);
  }

  // Mark released
  array->release = NULL;
}

static ArrowErrorCode ArrowArraySetStorageType(struct ArrowArray* array,
                                               enum ArrowType storage_type) {
  switch (storage_type) {
    case NANOARROW_TYPE_UNINITIALIZED:
    case NANOARROW_TYPE_NA:
      array->n_buffers = 0;
      break;

    case NANOARROW_TYPE_FIXED_SIZE_LIST:
    case NANOARROW_TYPE_STRUCT:
    case NANOARROW_TYPE_SPARSE_UNION:
      array->n_buffers = 1;
      break;

    case NANOARROW_TYPE_LIST:
    case NANOARROW_TYPE_LARGE_LIST:
    case NANOARROW_TYPE_MAP:
    case NANOARROW_TYPE_BOOL:
    case NANOARROW_TYPE_UINT8:
    case NANOARROW_TYPE_INT8:
    case NANOARROW_TYPE_UINT16:
    case NANOARROW_TYPE_INT16:
    case NANOARROW_TYPE_UINT32:
    case NANOARROW_TYPE_INT32:
    case NANOARROW_TYPE_UINT64:
    case NANOARROW_TYPE_INT64:
    case NANOARROW_TYPE_HALF_FLOAT:
    case NANOARROW_TYPE_FLOAT:
    case NANOARROW_TYPE_DOUBLE:
    case NANOARROW_TYPE_DECIMAL128:
    case NANOARROW_TYPE_DECIMAL256:
    case NANOARROW_TYPE_INTERVAL_MONTHS:
    case NANOARROW_TYPE_INTERVAL_DAY_TIME:
    case NANOARROW_TYPE_INTERVAL_MONTH_DAY_NANO:
    case NANOARROW_TYPE_FIXED_SIZE_BINARY:
    case NANOARROW_TYPE_DENSE_UNION:
      array->n_buffers = 2;
      break;

    case NANOARROW_TYPE_STRING:
    case NANOARROW_TYPE_LARGE_STRING:
    case NANOARROW_TYPE_BINARY:
    case NANOARROW_TYPE_LARGE_BINARY:
      array->n_buffers = 3;
      break;

    default:
      return EINVAL;

      return NANOARROW_OK;
  }

  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;
  private_data->storage_type = storage_type;
  return NANOARROW_OK;
}

ArrowErrorCode ArrowArrayInitFromType(struct ArrowArray* array,
                                      enum ArrowType storage_type) {
  array->length = 0;
  array->null_count = 0;
  array->offset = 0;
  array->n_buffers = 0;
  array->n_children = 0;
  array->buffers = NULL;
  array->children = NULL;
  array->dictionary = NULL;
  array->release = &ArrowArrayRelease;
  array->private_data = NULL;

  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)ArrowMalloc(sizeof(struct ArrowArrayPrivateData));
  if (private_data == NULL) {
    array->release = NULL;
    return ENOMEM;
  }

  ArrowBitmapInit(&private_data->bitmap);
  ArrowBufferInit(&private_data->buffers[0]);
  ArrowBufferInit(&private_data->buffers[1]);
  private_data->buffer_data[0] = NULL;
  private_data->buffer_data[1] = NULL;
  private_data->buffer_data[2] = NULL;

  array->private_data = private_data;
  array->buffers = (const void**)(&private_data->buffer_data);

  int result = ArrowArraySetStorageType(array, storage_type);
  if (result != NANOARROW_OK) {
    array->release(array);
    return result;
  }

  ArrowLayoutInit(&private_data->layout, storage_type);
  // We can only know this not to be true when initializing based on a schema
  // so assume this to be true.
  private_data->union_type_id_is_child_index = 1;
  return NANOARROW_OK;
}

ArrowErrorCode ArrowArrayInitFromArrayView(struct ArrowArray* array,
                                           struct ArrowArrayView* array_view,
                                           struct ArrowError* error) {
  NANOARROW_RETURN_NOT_OK_WITH_ERROR(
      ArrowArrayInitFromType(array, array_view->storage_type), error);
  int result;

  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;
  private_data->layout = array_view->layout;

  if (array_view->n_children > 0) {
    result = ArrowArrayAllocateChildren(array, array_view->n_children);
    if (result != NANOARROW_OK) {
      array->release(array);
      return result;
    }

    for (int64_t i = 0; i < array_view->n_children; i++) {
      result =
          ArrowArrayInitFromArrayView(array->children[i], array_view->children[i], error);
      if (result != NANOARROW_OK) {
        array->release(array);
        return result;
      }
    }
  }

  if (array_view->dictionary != NULL) {
    result = ArrowArrayAllocateDictionary(array);
    if (result != NANOARROW_OK) {
      array->release(array);
      return result;
    }

    result =
        ArrowArrayInitFromArrayView(array->dictionary, array_view->dictionary, error);
    if (result != NANOARROW_OK) {
      array->release(array);
      return result;
    }
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowArrayInitFromSchema(struct ArrowArray* array,
                                        struct ArrowSchema* schema,
                                        struct ArrowError* error) {
  struct ArrowArrayView array_view;
  NANOARROW_RETURN_NOT_OK(ArrowArrayViewInitFromSchema(&array_view, schema, error));
  NANOARROW_RETURN_NOT_OK(ArrowArrayInitFromArrayView(array, &array_view, error));
  if (array_view.storage_type == NANOARROW_TYPE_DENSE_UNION ||
      array_view.storage_type == NANOARROW_TYPE_SPARSE_UNION) {
    struct ArrowArrayPrivateData* private_data =
        (struct ArrowArrayPrivateData*)array->private_data;
    // We can still build arrays if this isn't true; however, the append
    // functions won't work. Instead, we store this value and error only
    // when StartAppending is called.
    private_data->union_type_id_is_child_index =
        _ArrowUnionTypeIdsWillEqualChildIndices(schema->format + 4, schema->n_children);
  }

  ArrowArrayViewReset(&array_view);
  return NANOARROW_OK;
}

ArrowErrorCode ArrowArrayAllocateChildren(struct ArrowArray* array, int64_t n_children) {
  if (array->children != NULL) {
    return EINVAL;
  }

  if (n_children == 0) {
    return NANOARROW_OK;
  }

  array->children =
      (struct ArrowArray**)ArrowMalloc(n_children * sizeof(struct ArrowArray*));
  if (array->children == NULL) {
    return ENOMEM;
  }

  memset(array->children, 0, n_children * sizeof(struct ArrowArray*));

  for (int64_t i = 0; i < n_children; i++) {
    array->children[i] = (struct ArrowArray*)ArrowMalloc(sizeof(struct ArrowArray));
    if (array->children[i] == NULL) {
      return ENOMEM;
    }
    array->children[i]->release = NULL;
  }

  array->n_children = n_children;
  return NANOARROW_OK;
}

ArrowErrorCode ArrowArrayAllocateDictionary(struct ArrowArray* array) {
  if (array->dictionary != NULL) {
    return EINVAL;
  }

  array->dictionary = (struct ArrowArray*)ArrowMalloc(sizeof(struct ArrowArray));
  if (array->dictionary == NULL) {
    return ENOMEM;
  }

  array->dictionary->release = NULL;
  return NANOARROW_OK;
}

void ArrowArraySetValidityBitmap(struct ArrowArray* array, struct ArrowBitmap* bitmap) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;
  ArrowBufferMove(&bitmap->buffer, &private_data->bitmap.buffer);
  private_data->bitmap.size_bits = bitmap->size_bits;
  bitmap->size_bits = 0;
  private_data->buffer_data[0] = private_data->bitmap.buffer.data;
  array->null_count = -1;
}

ArrowErrorCode ArrowArraySetBuffer(struct ArrowArray* array, int64_t i,
                                   struct ArrowBuffer* buffer) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  switch (i) {
    case 0:
      ArrowBufferMove(buffer, &private_data->bitmap.buffer);
      private_data->buffer_data[i] = private_data->bitmap.buffer.data;
      break;
    case 1:
    case 2:
      ArrowBufferMove(buffer, &private_data->buffers[i - 1]);
      private_data->buffer_data[i] = private_data->buffers[i - 1].data;
      break;
    default:
      return EINVAL;
  }

  return NANOARROW_OK;
}

static ArrowErrorCode ArrowArrayViewInitFromArray(struct ArrowArrayView* array_view,
                                                  struct ArrowArray* array) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  ArrowArrayViewInitFromType(array_view, private_data->storage_type);
  array_view->layout = private_data->layout;
  array_view->array = array;
  array_view->length = array->length;
  array_view->offset = array->offset;
  array_view->null_count = array->null_count;

  array_view->buffer_views[0].data.as_uint8 = private_data->bitmap.buffer.data;
  array_view->buffer_views[0].size_bytes = private_data->bitmap.buffer.size_bytes;
  array_view->buffer_views[1].data.as_uint8 = private_data->buffers[0].data;
  array_view->buffer_views[1].size_bytes = private_data->buffers[0].size_bytes;
  array_view->buffer_views[2].data.as_uint8 = private_data->buffers[1].data;
  array_view->buffer_views[2].size_bytes = private_data->buffers[1].size_bytes;

  int result = ArrowArrayViewAllocateChildren(array_view, array->n_children);
  if (result != NANOARROW_OK) {
    ArrowArrayViewReset(array_view);
    return result;
  }

  for (int64_t i = 0; i < array->n_children; i++) {
    result = ArrowArrayViewInitFromArray(array_view->children[i], array->children[i]);
    if (result != NANOARROW_OK) {
      ArrowArrayViewReset(array_view);
      return result;
    }
  }

  if (array->dictionary != NULL) {
    result = ArrowArrayViewAllocateDictionary(array_view);
    if (result != NANOARROW_OK) {
      ArrowArrayViewReset(array_view);
      return result;
    }

    result = ArrowArrayViewInitFromArray(array_view->dictionary, array->dictionary);
    if (result != NANOARROW_OK) {
      ArrowArrayViewReset(array_view);
      return result;
    }
  }

  return NANOARROW_OK;
}

static ArrowErrorCode ArrowArrayReserveInternal(struct ArrowArray* array,
                                                struct ArrowArrayView* array_view) {
  // Loop through buffers and reserve the extra space that we know about
  for (int64_t i = 0; i < array->n_buffers; i++) {
    // Don't reserve on a validity buffer that hasn't been allocated yet
    if (array_view->layout.buffer_type[i] == NANOARROW_BUFFER_TYPE_VALIDITY &&
        ArrowArrayBuffer(array, i)->data == NULL) {
      continue;
    }

    int64_t additional_size_bytes =
        array_view->buffer_views[i].size_bytes - ArrowArrayBuffer(array, i)->size_bytes;

    if (additional_size_bytes > 0) {
      NANOARROW_RETURN_NOT_OK(
          ArrowBufferReserve(ArrowArrayBuffer(array, i), additional_size_bytes));
    }
  }

  // Recursively reserve children
  for (int64_t i = 0; i < array->n_children; i++) {
    NANOARROW_RETURN_NOT_OK(
        ArrowArrayReserveInternal(array->children[i], array_view->children[i]));
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowArrayReserve(struct ArrowArray* array,
                                 int64_t additional_size_elements) {
  struct ArrowArrayView array_view;
  NANOARROW_RETURN_NOT_OK(ArrowArrayViewInitFromArray(&array_view, array));

  // Calculate theoretical buffer sizes (recursively)
  ArrowArrayViewSetLength(&array_view, array->length + additional_size_elements);

  // Walk the structure (recursively)
  int result = ArrowArrayReserveInternal(array, &array_view);
  ArrowArrayViewReset(&array_view);
  if (result != NANOARROW_OK) {
    return result;
  }

  return NANOARROW_OK;
}

static ArrowErrorCode ArrowArrayFinalizeBuffers(struct ArrowArray* array) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  // The only buffer finalizing this currently does is make sure the data
  // buffer for (Large)String|Binary is never NULL
  switch (private_data->storage_type) {
    case NANOARROW_TYPE_BINARY:
    case NANOARROW_TYPE_STRING:
    case NANOARROW_TYPE_LARGE_BINARY:
    case NANOARROW_TYPE_LARGE_STRING:
      if (ArrowArrayBuffer(array, 2)->data == NULL) {
        ArrowBufferAppendUInt8(ArrowArrayBuffer(array, 2), 0);
      }
      break;
    default:
      break;
  }

  for (int64_t i = 0; i < array->n_children; i++) {
    NANOARROW_RETURN_NOT_OK(ArrowArrayFinalizeBuffers(array->children[i]));
  }

  if (array->dictionary != NULL) {
    NANOARROW_RETURN_NOT_OK(ArrowArrayFinalizeBuffers(array->dictionary));
  }

  return NANOARROW_OK;
}

static void ArrowArrayFlushInternalPointers(struct ArrowArray* array) {
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;

  for (int64_t i = 0; i < 3; i++) {
    private_data->buffer_data[i] = ArrowArrayBuffer(array, i)->data;
  }

  for (int64_t i = 0; i < array->n_children; i++) {
    ArrowArrayFlushInternalPointers(array->children[i]);
  }

  if (array->dictionary != NULL) {
    ArrowArrayFlushInternalPointers(array->dictionary);
  }
}

ArrowErrorCode ArrowArrayFinishBuilding(struct ArrowArray* array,
                                        enum ArrowValidationLevel validation_level,
                                        struct ArrowError* error) {
  // Even if the data buffer is size zero, the pointer value needed to be non-null
  // in some implementations (at least one version of Arrow C++ at the time this
  // was added). Only do this fix if we can assume CPU data access.
  if (validation_level >= NANOARROW_VALIDATION_LEVEL_DEFAULT) {
    NANOARROW_RETURN_NOT_OK_WITH_ERROR(ArrowArrayFinalizeBuffers(array), error);
  }

  // Make sure the value we get with array->buffers[i] is set to the actual
  // pointer (which may have changed from the original due to reallocation)
  ArrowArrayFlushInternalPointers(array);

  if (validation_level == NANOARROW_VALIDATION_LEVEL_NONE) {
    return NANOARROW_OK;
  }

  // For validation, initialize an ArrowArrayView with our known buffer sizes
  struct ArrowArrayView array_view;
  NANOARROW_RETURN_NOT_OK_WITH_ERROR(ArrowArrayViewInitFromArray(&array_view, array),
                                     error);
  int result = ArrowArrayViewValidate(&array_view, validation_level, error);
  ArrowArrayViewReset(&array_view);
  return result;
}

ArrowErrorCode ArrowArrayFinishBuildingDefault(struct ArrowArray* array,
                                               struct ArrowError* error) {
  return ArrowArrayFinishBuilding(array, NANOARROW_VALIDATION_LEVEL_DEFAULT, error);
}

void ArrowArrayViewInitFromType(struct ArrowArrayView* array_view,
                                enum ArrowType storage_type) {
  memset(array_view, 0, sizeof(struct ArrowArrayView));
  array_view->storage_type = storage_type;
  ArrowLayoutInit(&array_view->layout, storage_type);
}

ArrowErrorCode ArrowArrayViewAllocateChildren(struct ArrowArrayView* array_view,
                                              int64_t n_children) {
  if (array_view->children != NULL) {
    return EINVAL;
  }

  array_view->children =
      (struct ArrowArrayView**)ArrowMalloc(n_children * sizeof(struct ArrowArrayView*));
  if (array_view->children == NULL) {
    return ENOMEM;
  }

  for (int64_t i = 0; i < n_children; i++) {
    array_view->children[i] = NULL;
  }

  array_view->n_children = n_children;

  for (int64_t i = 0; i < n_children; i++) {
    array_view->children[i] =
        (struct ArrowArrayView*)ArrowMalloc(sizeof(struct ArrowArrayView));
    if (array_view->children[i] == NULL) {
      return ENOMEM;
    }
    ArrowArrayViewInitFromType(array_view->children[i], NANOARROW_TYPE_UNINITIALIZED);
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowArrayViewAllocateDictionary(struct ArrowArrayView* array_view) {
  if (array_view->dictionary != NULL) {
    return EINVAL;
  }

  array_view->dictionary =
      (struct ArrowArrayView*)ArrowMalloc(sizeof(struct ArrowArrayView));
  if (array_view->dictionary == NULL) {
    return ENOMEM;
  }

  ArrowArrayViewInitFromType(array_view->dictionary, NANOARROW_TYPE_UNINITIALIZED);
  return NANOARROW_OK;
}

ArrowErrorCode ArrowArrayViewInitFromSchema(struct ArrowArrayView* array_view,
                                            struct ArrowSchema* schema,
                                            struct ArrowError* error) {
  struct ArrowSchemaView schema_view;
  int result = ArrowSchemaViewInit(&schema_view, schema, error);
  if (result != NANOARROW_OK) {
    return result;
  }

  ArrowArrayViewInitFromType(array_view, schema_view.storage_type);
  array_view->layout = schema_view.layout;

  result = ArrowArrayViewAllocateChildren(array_view, schema->n_children);
  if (result != NANOARROW_OK) {
    ArrowErrorSet(error, "ArrowArrayViewAllocateChildren() failed");
    ArrowArrayViewReset(array_view);
    return result;
  }

  for (int64_t i = 0; i < schema->n_children; i++) {
    result =
        ArrowArrayViewInitFromSchema(array_view->children[i], schema->children[i], error);
    if (result != NANOARROW_OK) {
      ArrowArrayViewReset(array_view);
      return result;
    }
  }

  if (schema->dictionary != NULL) {
    result = ArrowArrayViewAllocateDictionary(array_view);
    if (result != NANOARROW_OK) {
      ArrowArrayViewReset(array_view);
      return result;
    }

    result =
        ArrowArrayViewInitFromSchema(array_view->dictionary, schema->dictionary, error);
    if (result != NANOARROW_OK) {
      ArrowArrayViewReset(array_view);
      return result;
    }
  }

  if (array_view->storage_type == NANOARROW_TYPE_SPARSE_UNION ||
      array_view->storage_type == NANOARROW_TYPE_DENSE_UNION) {
    array_view->union_type_id_map = (int8_t*)ArrowMalloc(256 * sizeof(int8_t));
    if (array_view->union_type_id_map == NULL) {
      return ENOMEM;
    }

    memset(array_view->union_type_id_map, -1, 256);
    int8_t n_type_ids = _ArrowParseUnionTypeIds(schema_view.union_type_ids,
                                                array_view->union_type_id_map + 128);
    for (int8_t child_index = 0; child_index < n_type_ids; child_index++) {
      int8_t type_id = array_view->union_type_id_map[128 + child_index];
      array_view->union_type_id_map[type_id] = child_index;
    }
  }

  return NANOARROW_OK;
}

void ArrowArrayViewReset(struct ArrowArrayView* array_view) {
  if (array_view->children != NULL) {
    for (int64_t i = 0; i < array_view->n_children; i++) {
      if (array_view->children[i] != NULL) {
        ArrowArrayViewReset(array_view->children[i]);
        ArrowFree(array_view->children[i]);
      }
    }

    ArrowFree(array_view->children);
  }

  if (array_view->dictionary != NULL) {
    ArrowArrayViewReset(array_view->dictionary);
    ArrowFree(array_view->dictionary);
  }

  if (array_view->union_type_id_map != NULL) {
    ArrowFree(array_view->union_type_id_map);
  }

  ArrowArrayViewInitFromType(array_view, NANOARROW_TYPE_UNINITIALIZED);
}

void ArrowArrayViewSetLength(struct ArrowArrayView* array_view, int64_t length) {
  for (int i = 0; i < 3; i++) {
    int64_t element_size_bytes = array_view->layout.element_size_bits[i] / 8;

    switch (array_view->layout.buffer_type[i]) {
      case NANOARROW_BUFFER_TYPE_VALIDITY:
        array_view->buffer_views[i].size_bytes = _ArrowBytesForBits(length);
        continue;
      case NANOARROW_BUFFER_TYPE_DATA_OFFSET:
        // Probably don't want/need to rely on the producer to have allocated an
        // offsets buffer of length 1 for a zero-size array
        array_view->buffer_views[i].size_bytes =
            (length != 0) * element_size_bytes * (length + 1);
        continue;
      case NANOARROW_BUFFER_TYPE_DATA:
        array_view->buffer_views[i].size_bytes =
            _ArrowRoundUpToMultipleOf8(array_view->layout.element_size_bits[i] * length) /
            8;
        continue;
      case NANOARROW_BUFFER_TYPE_TYPE_ID:
      case NANOARROW_BUFFER_TYPE_UNION_OFFSET:
        array_view->buffer_views[i].size_bytes = element_size_bytes * length;
        continue;
      case NANOARROW_BUFFER_TYPE_NONE:
        array_view->buffer_views[i].size_bytes = 0;
        continue;
    }
  }

  switch (array_view->storage_type) {
    case NANOARROW_TYPE_STRUCT:
    case NANOARROW_TYPE_SPARSE_UNION:
      for (int64_t i = 0; i < array_view->n_children; i++) {
        ArrowArrayViewSetLength(array_view->children[i], length);
      }
      break;
    case NANOARROW_TYPE_FIXED_SIZE_LIST:
      if (array_view->n_children >= 1) {
        ArrowArrayViewSetLength(array_view->children[0],
                                length * array_view->layout.child_size_elements);
      }
    default:
      break;
  }
}

// This version recursively extracts information from the array and stores it
// in the array view, performing any checks that require the original array.
static int ArrowArrayViewSetArrayInternal(struct ArrowArrayView* array_view,
                                          struct ArrowArray* array,
                                          struct ArrowError* error) {
  // Check length and offset
  if (array->offset < 0) {
    ArrowErrorSet(error, "Expected array offset >= 0 but found array offset of %ld",
                  (long)array->offset);
    return EINVAL;
  }

  if (array->length < 0) {
    ArrowErrorSet(error, "Expected array length >= 0 but found array length of %ld",
                  (long)array->length);
    return EINVAL;
  }

  array_view->array = array;
  array_view->offset = array->offset;
  array_view->length = array->length;
  array_view->null_count = array->null_count;

  int64_t buffers_required = 0;
  for (int i = 0; i < 3; i++) {
    if (array_view->layout.buffer_type[i] == NANOARROW_BUFFER_TYPE_NONE) {
      break;
    }

    buffers_required++;

    // Set buffer pointer
    array_view->buffer_views[i].data.data = array->buffers[i];

    // If non-null, set buffer size to unknown.
    if (array->buffers[i] == NULL) {
      array_view->buffer_views[i].size_bytes = 0;
    } else {
      array_view->buffer_views[i].size_bytes = -1;
    }
  }

  // Check the number of buffers
  if (buffers_required != array->n_buffers) {
    ArrowErrorSet(error, "Expected array with %d buffer(s) but found %d buffer(s)",
                  (int)buffers_required, (int)array->n_buffers);
    return EINVAL;
  }

  // Check number of children
  if (array_view->n_children != array->n_children) {
    ArrowErrorSet(error, "Expected %ld children but found %ld children",
                  (long)array_view->n_children, (long)array->n_children);
    return EINVAL;
  }

  // Recurse for children
  for (int64_t i = 0; i < array_view->n_children; i++) {
    NANOARROW_RETURN_NOT_OK(ArrowArrayViewSetArrayInternal(array_view->children[i],
                                                           array->children[i], error));
  }

  // Check dictionary
  if (array->dictionary == NULL && array_view->dictionary != NULL) {
    ArrowErrorSet(error, "Expected dictionary but found NULL");
    return EINVAL;
  }

  if (array->dictionary != NULL && array_view->dictionary == NULL) {
    ArrowErrorSet(error, "Expected NULL dictionary but found dictionary member");
    return EINVAL;
  }

  if (array->dictionary != NULL) {
    NANOARROW_RETURN_NOT_OK(
        ArrowArrayViewSetArrayInternal(array_view->dictionary, array->dictionary, error));
  }

  return NANOARROW_OK;
}

static int ArrowArrayViewValidateMinimal(struct ArrowArrayView* array_view,
                                         struct ArrowError* error) {
  // Calculate buffer sizes that do not require buffer access. If marked as
  // unknown, assign the buffer size; otherwise, validate it.
  int64_t offset_plus_length = array_view->offset + array_view->length;

  // Only loop over the first two buffers because the size of the third buffer
  // is always data dependent for all current Arrow types.
  for (int i = 0; i < 2; i++) {
    int64_t element_size_bytes = array_view->layout.element_size_bits[i] / 8;
    // Initialize with a value that will cause an error if accidentally used uninitialized
    int64_t min_buffer_size_bytes = array_view->buffer_views[i].size_bytes + 1;

    switch (array_view->layout.buffer_type[i]) {
      case NANOARROW_BUFFER_TYPE_VALIDITY:
        if (array_view->null_count == 0 && array_view->buffer_views[i].size_bytes == 0) {
          continue;
        }

        min_buffer_size_bytes = _ArrowBytesForBits(offset_plus_length);
        break;
      case NANOARROW_BUFFER_TYPE_DATA_OFFSET:
        // Probably don't want/need to rely on the producer to have allocated an
        // offsets buffer of length 1 for a zero-size array
        min_buffer_size_bytes =
            (offset_plus_length != 0) * element_size_bytes * (offset_plus_length + 1);
        break;
      case NANOARROW_BUFFER_TYPE_DATA:
        min_buffer_size_bytes =
            _ArrowRoundUpToMultipleOf8(array_view->layout.element_size_bits[i] *
                                       offset_plus_length) /
            8;
        break;
      case NANOARROW_BUFFER_TYPE_TYPE_ID:
      case NANOARROW_BUFFER_TYPE_UNION_OFFSET:
        min_buffer_size_bytes = element_size_bytes * offset_plus_length;
        break;
      case NANOARROW_BUFFER_TYPE_NONE:
        continue;
    }

    // Assign or validate buffer size
    if (array_view->buffer_views[i].size_bytes == -1) {
      array_view->buffer_views[i].size_bytes = min_buffer_size_bytes;
    } else if (array_view->buffer_views[i].size_bytes < min_buffer_size_bytes) {
      ArrowErrorSet(error,
                    "Expected %s array buffer %d to have size >= %ld bytes but found "
                    "buffer with %ld bytes",
                    ArrowTypeString(array_view->storage_type), (int)i,
                    (long)min_buffer_size_bytes,
                    (long)array_view->buffer_views[i].size_bytes);
      return EINVAL;
    }
  }

  // For list, fixed-size list and map views, we can validate the number of children
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_LIST:
    case NANOARROW_TYPE_LARGE_LIST:
    case NANOARROW_TYPE_FIXED_SIZE_LIST:
    case NANOARROW_TYPE_MAP:
      if (array_view->n_children != 1) {
        ArrowErrorSet(error, "Expected 1 child of %s array but found %ld child arrays",
                      ArrowTypeString(array_view->storage_type),
                      (long)array_view->n_children);
        return EINVAL;
      }
    default:
      break;
  }

  // For struct, the sparse union, and the fixed-size list views, we can validate child
  // lengths.
  int64_t child_min_length;
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_SPARSE_UNION:
    case NANOARROW_TYPE_STRUCT:
      child_min_length = (array_view->offset + array_view->length);
      for (int64_t i = 0; i < array_view->n_children; i++) {
        if (array_view->children[i]->length < child_min_length) {
          ArrowErrorSet(
              error,
              "Expected struct child %d to have length >= %ld but found child with "
              "length %ld",
              (int)(i + 1), (long)(child_min_length),
              (long)array_view->children[i]->length);
          return EINVAL;
        }
      }
      break;

    case NANOARROW_TYPE_FIXED_SIZE_LIST:
      child_min_length = (array_view->offset + array_view->length) *
                         array_view->layout.child_size_elements;
      if (array_view->children[0]->length < child_min_length) {
        ArrowErrorSet(error,
                      "Expected child of fixed_size_list array to have length >= %ld but "
                      "found array with length %ld",
                      (long)child_min_length, (long)array_view->children[0]->length);
        return EINVAL;
      }
      break;
    default:
      break;
  }

  // Recurse for children
  for (int64_t i = 0; i < array_view->n_children; i++) {
    NANOARROW_RETURN_NOT_OK(
        ArrowArrayViewValidateMinimal(array_view->children[i], error));
  }

  // Recurse for dictionary
  if (array_view->dictionary != NULL) {
    NANOARROW_RETURN_NOT_OK(ArrowArrayViewValidateMinimal(array_view->dictionary, error));
  }

  return NANOARROW_OK;
}

static int ArrowArrayViewValidateDefault(struct ArrowArrayView* array_view,
                                         struct ArrowError* error) {
  // Perform minimal validation. This will validate or assign
  // buffer sizes as long as buffer access is not required.
  NANOARROW_RETURN_NOT_OK(ArrowArrayViewValidateMinimal(array_view, error));

  // Calculate buffer sizes or child lengths that require accessing the offsets
  // buffer. Where appropriate, validate that the first offset is >= 0.
  // If a buffer size is marked as unknown, assign it; otherwise, validate it.
  int64_t offset_plus_length = array_view->offset + array_view->length;

  int64_t first_offset;
  int64_t last_offset;
  switch (array_view->storage_type) {
    case NANOARROW_TYPE_STRING:
    case NANOARROW_TYPE_BINARY:
      if (array_view->buffer_views[1].size_bytes != 0) {
        first_offset = array_view->buffer_views[1].data.as_int32[0];
        if (first_offset < 0) {
          ArrowErrorSet(error, "Expected first offset >= 0 but found %ld",
                        (long)first_offset);
          return EINVAL;
        }

        last_offset = array_view->buffer_views[1].data.as_int32[offset_plus_length];

        // If the data buffer size is unknown, assign it; otherwise, check it
        if (array_view->buffer_views[2].size_bytes == -1) {
          array_view->buffer_views[2].size_bytes = last_offset;
        } else if (array_view->buffer_views[2].size_bytes < last_offset) {
          ArrowErrorSet(error,
                        "Expected %s array buffer 2 to have size >= %ld bytes but found "
                        "buffer with %ld bytes",
                        ArrowTypeString(array_view->storage_type), (long)last_offset,
                        (long)array_view->buffer_views[2].size_bytes);
          return EINVAL;
        }
      }
      break;

    case NANOARROW_TYPE_LARGE_STRING:
    case NANOARROW_TYPE_LARGE_BINARY:
      if (array_view->buffer_views[1].size_bytes != 0) {
        first_offset = array_view->buffer_views[1].data.as_int64[0];
        if (first_offset < 0) {
          ArrowErrorSet(error, "Expected first offset >= 0 but found %ld",
                        (long)first_offset);
          return EINVAL;
        }

        last_offset = array_view->buffer_views[1].data.as_int64[offset_plus_length];

        // If the data buffer size is unknown, assign it; otherwise, check it
        if (array_view->buffer_views[2].size_bytes == -1) {
          array_view->buffer_views[2].size_bytes = last_offset;
        } else if (array_view->buffer_views[2].size_bytes < last_offset) {
          ArrowErrorSet(error,
                        "Expected %s array buffer 2 to have size >= %ld bytes but found "
                        "buffer with %ld bytes",
                        ArrowTypeString(array_view->storage_type), (long)last_offset,
                        (long)array_view->buffer_views[2].size_bytes);
          return EINVAL;
        }
      }
      break;

    case NANOARROW_TYPE_STRUCT:
      for (int64_t i = 0; i < array_view->n_children; i++) {
        if (array_view->children[i]->length < offset_plus_length) {
          ArrowErrorSet(
              error,
              "Expected struct child %d to have length >= %ld but found child with "
              "length %ld",
              (int)(i + 1), (long)offset_plus_length,
              (long)array_view->children[i]->length);
          return EINVAL;
        }
      }
      break;

    case NANOARROW_TYPE_LIST:
    case NANOARROW_TYPE_MAP:
      if (array_view->buffer_views[1].size_bytes != 0) {
        first_offset = array_view->buffer_views[1].data.as_int32[0];
        if (first_offset < 0) {
          ArrowErrorSet(error, "Expected first offset >= 0 but found %ld",
                        (long)first_offset);
          return EINVAL;
        }

        last_offset = array_view->buffer_views[1].data.as_int32[offset_plus_length];
        if (array_view->children[0]->length < last_offset) {
          ArrowErrorSet(
              error,
              "Expected child of %s array to have length >= %ld but found array with "
              "length %ld",
              ArrowTypeString(array_view->storage_type), (long)last_offset,
              (long)array_view->children[0]->length);
          return EINVAL;
        }
      }
      break;

    case NANOARROW_TYPE_LARGE_LIST:
      if (array_view->buffer_views[1].size_bytes != 0) {
        first_offset = array_view->buffer_views[1].data.as_int64[0];
        if (first_offset < 0) {
          ArrowErrorSet(error, "Expected first offset >= 0 but found %ld",
                        (long)first_offset);
          return EINVAL;
        }

        last_offset = array_view->buffer_views[1].data.as_int64[offset_plus_length];
        if (array_view->children[0]->length < last_offset) {
          ArrowErrorSet(
              error,
              "Expected child of large list array to have length >= %ld but found array "
              "with length %ld",
              (long)last_offset, (long)array_view->children[0]->length);
          return EINVAL;
        }
      }
      break;
    default:
      break;
  }

  // Recurse for children
  for (int64_t i = 0; i < array_view->n_children; i++) {
    NANOARROW_RETURN_NOT_OK(
        ArrowArrayViewValidateDefault(array_view->children[i], error));
  }

  // Recurse for dictionary
  if (array_view->dictionary != NULL) {
    NANOARROW_RETURN_NOT_OK(ArrowArrayViewValidateDefault(array_view->dictionary, error));
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowArrayViewSetArray(struct ArrowArrayView* array_view,
                                      struct ArrowArray* array,
                                      struct ArrowError* error) {
  // Extract information from the array into the array view
  NANOARROW_RETURN_NOT_OK(ArrowArrayViewSetArrayInternal(array_view, array, error));

  // Run default validation. Because we've marked all non-NULL buffers as having unknown
  // size, validation will also update the buffer sizes as it goes.
  NANOARROW_RETURN_NOT_OK(ArrowArrayViewValidateDefault(array_view, error));

  return NANOARROW_OK;
}

ArrowErrorCode ArrowArrayViewSetArrayMinimal(struct ArrowArrayView* array_view,
                                             struct ArrowArray* array,
                                             struct ArrowError* error) {
  // Extract information from the array into the array view
  NANOARROW_RETURN_NOT_OK(ArrowArrayViewSetArrayInternal(array_view, array, error));

  // Run default validation. Because we've marked all non-NULL buffers as having unknown
  // size, validation will also update the buffer sizes as it goes.
  NANOARROW_RETURN_NOT_OK(ArrowArrayViewValidateMinimal(array_view, error));

  return NANOARROW_OK;
}

static int ArrowAssertIncreasingInt32(struct ArrowBufferView view,
                                      struct ArrowError* error) {
  if (view.size_bytes <= (int64_t)sizeof(int32_t)) {
    return NANOARROW_OK;
  }

  for (int64_t i = 1; i < view.size_bytes / (int64_t)sizeof(int32_t); i++) {
    int32_t diff = view.data.as_int32[i] - view.data.as_int32[i - 1];
    if (diff < 0) {
      ArrowErrorSet(error, "[%ld] Expected element size >= 0 but found element size %ld",
                    (long)i, (long)diff);
      return EINVAL;
    }
  }

  return NANOARROW_OK;
}

static int ArrowAssertIncreasingInt64(struct ArrowBufferView view,
                                      struct ArrowError* error) {
  if (view.size_bytes <= (int64_t)sizeof(int64_t)) {
    return NANOARROW_OK;
  }

  for (int64_t i = 1; i < view.size_bytes / (int64_t)sizeof(int64_t); i++) {
    int64_t diff = view.data.as_int64[i] - view.data.as_int64[i - 1];
    if (diff < 0) {
      ArrowErrorSet(error, "[%ld] Expected element size >= 0 but found element size %ld",
                    (long)i, (long)diff);
      return EINVAL;
    }
  }

  return NANOARROW_OK;
}

static int ArrowAssertRangeInt8(struct ArrowBufferView view, int8_t min_value,
                                int8_t max_value, struct ArrowError* error) {
  for (int64_t i = 0; i < view.size_bytes; i++) {
    if (view.data.as_int8[i] < min_value || view.data.as_int8[i] > max_value) {
      ArrowErrorSet(error,
                    "[%ld] Expected buffer value between %d and %d but found value %d",
                    (long)i, (int)min_value, (int)max_value, (int)view.data.as_int8[i]);
      return EINVAL;
    }
  }

  return NANOARROW_OK;
}

static int ArrowAssertInt8In(struct ArrowBufferView view, const int8_t* values,
                             int64_t n_values, struct ArrowError* error) {
  for (int64_t i = 0; i < view.size_bytes; i++) {
    int item_found = 0;
    for (int64_t j = 0; j < n_values; j++) {
      if (view.data.as_int8[i] == values[j]) {
        item_found = 1;
        break;
      }
    }

    if (!item_found) {
      ArrowErrorSet(error, "[%ld] Unexpected buffer value %d", (long)i,
                    (int)view.data.as_int8[i]);
      return EINVAL;
    }
  }

  return NANOARROW_OK;
}

static int ArrowArrayViewValidateFull(struct ArrowArrayView* array_view,
                                      struct ArrowError* error) {
  for (int i = 0; i < 3; i++) {
    switch (array_view->layout.buffer_type[i]) {
      case NANOARROW_BUFFER_TYPE_DATA_OFFSET:
        if (array_view->layout.element_size_bits[i] == 32) {
          NANOARROW_RETURN_NOT_OK(
              ArrowAssertIncreasingInt32(array_view->buffer_views[i], error));
        } else {
          NANOARROW_RETURN_NOT_OK(
              ArrowAssertIncreasingInt64(array_view->buffer_views[i], error));
        }
        break;
      default:
        break;
    }
  }

  if (array_view->storage_type == NANOARROW_TYPE_DENSE_UNION ||
      array_view->storage_type == NANOARROW_TYPE_SPARSE_UNION) {
    if (array_view->union_type_id_map == NULL) {
      // If the union_type_id map is NULL (e.g., when using ArrowArrayInitFromType() +
      // ArrowArrayAllocateChildren() + ArrowArrayFinishBuilding()), we don't have enough
      // information to validate this buffer.
      ArrowErrorSet(error,
                    "Insufficient information provided for validation of union array");
      return EINVAL;
    } else if (_ArrowParsedUnionTypeIdsWillEqualChildIndices(
                   array_view->union_type_id_map, array_view->n_children,
                   array_view->n_children)) {
      NANOARROW_RETURN_NOT_OK(ArrowAssertRangeInt8(
          array_view->buffer_views[0], 0, (int8_t)(array_view->n_children - 1), error));
    } else {
      NANOARROW_RETURN_NOT_OK(ArrowAssertInt8In(array_view->buffer_views[0],
                                                array_view->union_type_id_map + 128,
                                                array_view->n_children, error));
    }
  }

  if (array_view->storage_type == NANOARROW_TYPE_DENSE_UNION &&
      array_view->union_type_id_map != NULL) {
    // Check that offsets refer to child elements that actually exist
    for (int64_t i = 0; i < array_view->length; i++) {
      int8_t child_id = ArrowArrayViewUnionChildIndex(array_view, i);
      int64_t offset = ArrowArrayViewUnionChildOffset(array_view, i);
      int64_t child_length = array_view->children[child_id]->length;
      if (offset < 0 || offset > child_length) {
        ArrowErrorSet(
            error,
            "[%ld] Expected union offset for child id %d to be between 0 and %ld but "
            "found offset value %ld",
            (long)i, (int)child_id, (long)child_length, offset);
        return EINVAL;
      }
    }
  }

  // Recurse for children
  for (int64_t i = 0; i < array_view->n_children; i++) {
    NANOARROW_RETURN_NOT_OK(ArrowArrayViewValidateFull(array_view->children[i], error));
  }

  // Dictionary valiation not implemented
  if (array_view->dictionary != NULL) {
    ArrowErrorSet(error, "Validation for dictionary-encoded arrays is not implemented");
    return ENOTSUP;
  }

  return NANOARROW_OK;
}

ArrowErrorCode ArrowArrayViewValidate(struct ArrowArrayView* array_view,
                                      enum ArrowValidationLevel validation_level,
                                      struct ArrowError* error) {
  switch (validation_level) {
    case NANOARROW_VALIDATION_LEVEL_NONE:
      return NANOARROW_OK;
    case NANOARROW_VALIDATION_LEVEL_MINIMAL:
      return ArrowArrayViewValidateMinimal(array_view, error);
    case NANOARROW_VALIDATION_LEVEL_DEFAULT:
      return ArrowArrayViewValidateDefault(array_view, error);
    case NANOARROW_VALIDATION_LEVEL_FULL:
      NANOARROW_RETURN_NOT_OK(ArrowArrayViewValidateDefault(array_view, error));
      return ArrowArrayViewValidateFull(array_view, error);
  }

  ArrowErrorSet(error, "validation_level not recognized");
  return EINVAL;
}
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

#include <errno.h>

#include "nanoarrow.h"

struct BasicArrayStreamPrivate {
  struct ArrowSchema schema;
  int64_t n_arrays;
  struct ArrowArray* arrays;
  int64_t arrays_i;
};

static int ArrowBasicArrayStreamGetSchema(struct ArrowArrayStream* array_stream,
                                          struct ArrowSchema* schema) {
  if (array_stream == NULL || array_stream->release == NULL) {
    return EINVAL;
  }

  struct BasicArrayStreamPrivate* private_data =
      (struct BasicArrayStreamPrivate*)array_stream->private_data;
  return ArrowSchemaDeepCopy(&private_data->schema, schema);
}

static int ArrowBasicArrayStreamGetNext(struct ArrowArrayStream* array_stream,
                                        struct ArrowArray* array) {
  if (array_stream == NULL || array_stream->release == NULL) {
    return EINVAL;
  }

  struct BasicArrayStreamPrivate* private_data =
      (struct BasicArrayStreamPrivate*)array_stream->private_data;

  if (private_data->arrays_i == private_data->n_arrays) {
    array->release = NULL;
    return NANOARROW_OK;
  }

  ArrowArrayMove(&private_data->arrays[private_data->arrays_i++], array);
  return NANOARROW_OK;
}

static const char* ArrowBasicArrayStreamGetLastError(
    struct ArrowArrayStream* array_stream) {
  return NULL;
}

static void ArrowBasicArrayStreamRelease(struct ArrowArrayStream* array_stream) {
  if (array_stream == NULL || array_stream->release == NULL) {
    return;
  }

  struct BasicArrayStreamPrivate* private_data =
      (struct BasicArrayStreamPrivate*)array_stream->private_data;

  if (private_data->schema.release != NULL) {
    private_data->schema.release(&private_data->schema);
  }

  for (int64_t i = 0; i < private_data->n_arrays; i++) {
    if (private_data->arrays[i].release != NULL) {
      private_data->arrays[i].release(&private_data->arrays[i]);
    }
  }

  if (private_data->arrays != NULL) {
    ArrowFree(private_data->arrays);
  }

  ArrowFree(private_data);
  array_stream->release = NULL;
}

ArrowErrorCode ArrowBasicArrayStreamInit(struct ArrowArrayStream* array_stream,
                                         struct ArrowSchema* schema, int64_t n_arrays) {
  struct BasicArrayStreamPrivate* private_data =
      (struct BasicArrayStreamPrivate*)ArrowMalloc(
          sizeof(struct BasicArrayStreamPrivate));
  if (private_data == NULL) {
    return ENOMEM;
  }

  ArrowSchemaMove(schema, &private_data->schema);

  private_data->n_arrays = n_arrays;
  private_data->arrays = NULL;
  private_data->arrays_i = 0;

  if (n_arrays > 0) {
    private_data->arrays =
        (struct ArrowArray*)ArrowMalloc(n_arrays * sizeof(struct ArrowArray));
    if (private_data->arrays == NULL) {
      ArrowBasicArrayStreamRelease(array_stream);
      return ENOMEM;
    }
  }

  for (int64_t i = 0; i < private_data->n_arrays; i++) {
    private_data->arrays[i].release = NULL;
  }

  array_stream->get_schema = &ArrowBasicArrayStreamGetSchema;
  array_stream->get_next = &ArrowBasicArrayStreamGetNext;
  array_stream->get_last_error = ArrowBasicArrayStreamGetLastError;
  array_stream->release = ArrowBasicArrayStreamRelease;
  array_stream->private_data = private_data;
  return NANOARROW_OK;
}

void ArrowBasicArrayStreamSetArray(struct ArrowArrayStream* array_stream, int64_t i,
                                   struct ArrowArray* array) {
  struct BasicArrayStreamPrivate* private_data =
      (struct BasicArrayStreamPrivate*)array_stream->private_data;
  ArrowArrayMove(array, &private_data->arrays[i]);
}

ArrowErrorCode ArrowBasicArrayStreamValidate(struct ArrowArrayStream* array_stream,
                                             struct ArrowError* error) {
  struct BasicArrayStreamPrivate* private_data =
      (struct BasicArrayStreamPrivate*)array_stream->private_data;

  struct ArrowArrayView array_view;
  NANOARROW_RETURN_NOT_OK(
      ArrowArrayViewInitFromSchema(&array_view, &private_data->schema, error));

  for (int64_t i = 0; i < private_data->n_arrays; i++) {
    if (private_data->arrays[i].release != NULL) {
      int result = ArrowArrayViewSetArray(&array_view, &private_data->arrays[i], error);
      if (result != NANOARROW_OK) {
        ArrowArrayViewReset(&array_view);
        return result;
      }
    }
  }

  ArrowArrayViewReset(&array_view);
  return NANOARROW_OK;
}
