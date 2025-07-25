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

/// \file abi.h Arrow C Data Interface
///
/// The Arrow C Data interface defines a very small, stable set
/// of C definitions which can be easily copied into any project's
/// source code and vendored to be used for columnar data interchange
/// in the Arrow format. For non-C/C++ languages and runtimes,
/// it should be almost as easy to translate the C definitions into
/// the corresponding C FFI declarations.
///
/// Applications and libraries can therefore work with Arrow memory
/// without necessarily using the Arrow libraries or reinventing
/// the wheel. Developers can choose between tight integration
/// with the Arrow software project or minimal integration with
/// the Arrow format only.

#pragma once

#include <stdint.h>

// Spec and documentation: https://arrow.apache.org/docs/format/CDataInterface.html

#ifdef __cplusplus
extern "C" {
#endif

#ifndef ARROW_C_DATA_INTERFACE
#  define ARROW_C_DATA_INTERFACE

#  define ARROW_FLAG_DICTIONARY_ORDERED 1
#  define ARROW_FLAG_NULLABLE 2
#  define ARROW_FLAG_MAP_KEYS_SORTED 4

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

#  define ARROW_STATISTICS_KEY_AVERAGE_BYTE_WIDTH_EXACT "ARROW:average_byte_width:exact"
#  define ARROW_STATISTICS_KEY_AVERAGE_BYTE_WIDTH_APPROXIMATE \
    "ARROW:average_byte_width:approximate"
#  define ARROW_STATISTICS_KEY_DISTINCT_COUNT_EXACT "ARROW:distinct_count:exact"
#  define ARROW_STATISTICS_KEY_DISTINCT_COUNT_APPROXIMATE \
    "ARROW:distinct_count:approximate"
#  define ARROW_STATISTICS_KEY_MAX_BYTE_WIDTH_EXACT "ARROW:max_byte_width:exact"
#  define ARROW_STATISTICS_KEY_MAX_BYTE_WIDTH_APPROXIMATE \
    "ARROW:max_byte_width:approximate"
#  define ARROW_STATISTICS_KEY_MAX_VALUE_EXACT "ARROW:max_value:exact"
#  define ARROW_STATISTICS_KEY_MAX_VALUE_APPROXIMATE "ARROW:max_value:approximate"
#  define ARROW_STATISTICS_KEY_MIN_VALUE_EXACT "ARROW:min_value:exact"
#  define ARROW_STATISTICS_KEY_MIN_VALUE_APPROXIMATE "ARROW:min_value:approximate"
#  define ARROW_STATISTICS_KEY_NULL_COUNT_EXACT "ARROW:null_count:exact"
#  define ARROW_STATISTICS_KEY_NULL_COUNT_APPROXIMATE "ARROW:null_count:approximate"
#  define ARROW_STATISTICS_KEY_ROW_COUNT_EXACT "ARROW:row_count:exact"
#  define ARROW_STATISTICS_KEY_ROW_COUNT_APPROXIMATE "ARROW:row_count:approximate"

#endif  // ARROW_C_DATA_INTERFACE

#ifndef ARROW_C_DEVICE_DATA_INTERFACE
#  define ARROW_C_DEVICE_DATA_INTERFACE

// Spec and Documentation: https://arrow.apache.org/docs/format/CDeviceDataInterface.html

// DeviceType for the allocated memory
typedef int32_t ArrowDeviceType;

// CPU device, same as using ArrowArray directly
#  define ARROW_DEVICE_CPU 1
// CUDA GPU Device
#  define ARROW_DEVICE_CUDA 2
// Pinned CUDA CPU memory by cudaMallocHost
#  define ARROW_DEVICE_CUDA_HOST 3
// OpenCL Device
#  define ARROW_DEVICE_OPENCL 4
// Vulkan buffer for next-gen graphics
#  define ARROW_DEVICE_VULKAN 7
// Metal for Apple GPU
#  define ARROW_DEVICE_METAL 8
// Verilog simulator buffer
#  define ARROW_DEVICE_VPI 9
// ROCm GPUs for AMD GPUs
#  define ARROW_DEVICE_ROCM 10
// Pinned ROCm CPU memory allocated by hipMallocHost
#  define ARROW_DEVICE_ROCM_HOST 11
// Reserved for extension
#  define ARROW_DEVICE_EXT_DEV 12
// CUDA managed/unified memory allocated by cudaMallocManaged
#  define ARROW_DEVICE_CUDA_MANAGED 13
// unified shared memory allocated on a oneAPI non-partitioned device.
#  define ARROW_DEVICE_ONEAPI 14
// GPU support for next-gen WebGPU standard
#  define ARROW_DEVICE_WEBGPU 15
// Qualcomm Hexagon DSP
#  define ARROW_DEVICE_HEXAGON 16

struct ArrowDeviceArray {
  // the Allocated Array
  //
  // the buffers in the array (along with the buffers of any
  // children) are what is allocated on the device.
  struct ArrowArray array;
  // The device id to identify a specific device
  int64_t device_id;
  // The type of device which can access this memory.
  ArrowDeviceType device_type;
  // An event-like object to synchronize on if needed.
  void* sync_event;
  // Reserved bytes for future expansion.
  int64_t reserved[3];
};

#endif  // ARROW_C_DEVICE_DATA_INTERFACE

#ifndef ARROW_C_STREAM_INTERFACE
#  define ARROW_C_STREAM_INTERFACE

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

#ifndef ARROW_C_DEVICE_STREAM_INTERFACE
#  define ARROW_C_DEVICE_STREAM_INTERFACE

// Equivalent to ArrowArrayStream, but for ArrowDeviceArrays.
//
// This stream is intended to provide a stream of data on a single
// device, if a producer wants data to be produced on multiple devices
// then multiple streams should be provided. One per device.
struct ArrowDeviceArrayStream {
  // The device that this stream produces data on.
  ArrowDeviceType device_type;

  // Callback to get the stream schema
  // (will be the same for all arrays in the stream).
  //
  // Return value 0 if successful, an `errno`-compatible error code otherwise.
  //
  // If successful, the ArrowSchema must be released independently from the stream.
  // The schema should be accessible via CPU memory.
  int (*get_schema)(struct ArrowDeviceArrayStream* self, struct ArrowSchema* out);

  // Callback to get the next array
  // (if no error and the array is released, the stream has ended)
  //
  // Return value: 0 if successful, an `errno`-compatible error code otherwise.
  //
  // If successful, the ArrowDeviceArray must be released independently from the stream.
  int (*get_next)(struct ArrowDeviceArrayStream* self, struct ArrowDeviceArray* out);

  // Callback to get optional detailed error information.
  // This must only be called if the last stream operation failed
  // with a non-0 return code.
  //
  // Return value: pointer to a null-terminated character array describing
  // the last error, or NULL if no description is available.
  //
  // The returned pointer is only valid until the next operation on this stream
  // (including release).
  const char* (*get_last_error)(struct ArrowDeviceArrayStream* self);

  // Release callback: release the stream's own resources.
  // Note that arrays returned by `get_next` must be individually released.
  void (*release)(struct ArrowDeviceArrayStream* self);

  // Opaque producer-specific data
  void* private_data;
};

#endif  // ARROW_C_DEVICE_STREAM_INTERFACE

#ifndef ARROW_C_ASYNC_STREAM_INTERFACE
#  define ARROW_C_ASYNC_STREAM_INTERFACE

// EXPERIMENTAL: ArrowAsyncTask represents available data from a producer that was passed
// to an invocation of `on_next_task` on the ArrowAsyncDeviceStreamHandler.
//
// The reason for this Task approach instead of the Async interface returning
// the Array directly is to allow for more complex thread handling and reducing
// context switching and data transfers between CPU cores (e.g. from one L1/L2
// cache to another) if desired.
//
// For example, the `on_next_task` callback can be called when data is ready, while
// the producer puts potential "decoding" logic in the `ArrowAsyncTask` object. This
// allows for the producer to manage the I/O on one thread which calls `on_next_task`
// and the consumer can determine when the decoding (producer logic in the `extract_data`
// callback of the task) occurs and on which thread, to avoid a CPU core transfer
// (data staying in the L2 cache).
struct ArrowAsyncTask {
  // This callback should populate the ArrowDeviceArray associated with this task.
  // The order of ArrowAsyncTasks provided by the producer enables a consumer to
  // ensure the order of data to process.
  //
  // This function is expected to be synchronous, but should not perform any blocking
  // I/O. Ideally it should be as cheap as possible so as to not tie up the consumer
  // thread unnecessarily.
  //
  // Returns: 0 if successful, errno-compatible error otherwise.
  //
  // If a non-0 value is returned then it should be followed by a call to `on_error`
  // on the appropriate ArrowAsyncDeviceStreamHandler. This is because it's highly
  // likely that whatever is calling this function may be entirely disconnected from
  // the current control flow. Indicating an error here with a non-zero return allows
  // the current flow to be aware of the error occurring, while still allowing any
  // logging or error handling to still be centralized in the `on_error` callback of
  // the original Async handler.
  //
  // Rather than a release callback, any required cleanup should be performed as part
  // of the invocation of `extract_data`. Ownership of the Array is passed to the consumer
  // calling this, and so it must be released separately.
  //
  // It is only valid to call this method exactly once.
  int (*extract_data)(struct ArrowAsyncTask* self, struct ArrowDeviceArray* out);

  // opaque task-specific data
  void* private_data;
};

// EXPERIMENTAL: ArrowAsyncProducer represents a 1-to-1 relationship between an async
// producer and consumer. This object allows the consumer to perform backpressure and flow
// control on the asynchronous stream processing. This object must be owned by the
// producer who creates it, and thus is responsible for cleaning it up.
struct ArrowAsyncProducer {
  // The device type that this stream produces data on.
  ArrowDeviceType device_type;

  // A consumer must call this function to start receiving on_next_task calls.
  //
  // It *must* be valid to call this synchronously from within `on_next_task` or
  // `on_schema`, but this function *must not* immediately call `on_next_task` so as
  // to avoid recursion and reentrant callbacks.
  //
  // After cancel has been called, additional calls to this function must be NOPs,
  // but allowed. While not cancelled, calling this function must register the
  // given number of additional arrays/batches to be produced with the producer.
  // The producer should only call `on_next_task` at most the registered number
  // of arrays before propagating backpressure.
  //
  // Any error encountered by calling request must be propagated by calling the `on_error`
  // callback of the ArrowAsyncDeviceStreamHandler.
  //
  // While not cancelled, any subsequent calls to `on_next_task`, `on_error` or
  // `release` should be scheduled by the producer to be called later.
  //
  // It is invalid for a consumer to call this with a value of n <= 0, producers should
  // error if given such a value.
  void (*request)(struct ArrowAsyncProducer* self, int64_t n);

  // This cancel callback signals a producer that it must eventually stop making calls
  // to on_next_task. It must be idempotent and thread-safe. After calling cancel once,
  // subsequent calls must be NOPs. This must not call any consumer-side handlers other
  // than `on_error`.
  //
  // It is not required that calling cancel affect the producer immediately, only that it
  // must eventually stop calling on_next_task and subsequently call release on the
  // async handler. As such, a consumer must be prepared to receive one or more calls to
  // `on_next_task` even after calling cancel if there are still requested arrays pending.
  //
  // Successful cancellation should *not* result in the producer calling `on_error`, it
  // should finish out any remaining tasks and eventually call `release`.
  //
  // Any error encountered during handling a call to cancel must be reported via the
  // on_error callback on the async stream handler.
  void (*cancel)(struct ArrowAsyncProducer* self);

  // Any additional metadata tied to a specific stream of data. This must either be NULL
  // or a valid pointer to metadata which is encoded in the same way schema metadata
  // would be. Non-null metadata must be valid for the lifetime of this object. As an
  // example a producer could use this to provide the total number of rows and/or batches
  // in the stream if known.
  const char* additional_metadata;

  // producer-specific opaque data.
  void* private_data;
};

// EXPERIMENTAL: Similar to ArrowDeviceArrayStream, except designed for an asynchronous
// style of interaction. While ArrowDeviceArrayStream provides producer
// defined callbacks, this is intended to be created by the consumer instead.
// The consumer passes this handler to the producer, which in turn uses the
// callbacks to inform the consumer of events in the stream.
struct ArrowAsyncDeviceStreamHandler {
  // Handler for receiving a schema. The passed in stream_schema must be
  // released or moved by the handler (producer is giving ownership of the schema to
  // the handler, but not ownership of the top level object itself).
  //
  // With the exception of an error occurring (on_error), this must be the first
  // callback function which is called by a producer and must only be called exactly
  // once. As such, the producer should provide a valid ArrowAsyncProducer instance
  // so the consumer can control the flow. See the documentation on ArrowAsyncProducer
  // for how it works. The ArrowAsyncProducer is owned by the producer who calls this
  // function and thus the producer is responsible for cleaning it up when calling
  // the release callback of this handler.
  //
  // If there is any additional metadata tied to this stream, it will be provided as
  // a non-null value for the `additional_metadata` field of the ArrowAsyncProducer
  // which will be valid at least until the release callback is called.
  //
  // Return value: 0 if successful, `errno`-compatible error otherwise
  //
  // A producer that receives a non-zero return here should stop producing and eventually
  // call release instead.
  int (*on_schema)(struct ArrowAsyncDeviceStreamHandler* self,
                   struct ArrowSchema* stream_schema);

  // Handler for receiving data. This is called when data is available providing an
  // ArrowAsyncTask struct to signify it. The producer indicates the end of the stream
  // by passing NULL as the value for the task rather than a valid pointer to a task.
  // The task object is only valid for the lifetime of this function call, if a consumer
  // wants to utilize it after this function returns, it must copy or move the contents
  // of it to a new ArrowAsyncTask object.
  //
  // The `request` callback of a provided ArrowAsyncProducer must be called in order
  // to start receiving calls to this handler.
  //
  // The metadata argument can be null or can be used by a producer
  // to pass arbitrary extra information to the consumer (such as total number
  // of rows, context info, or otherwise). The data should be passed using the same
  // encoding as the metadata within the ArrowSchema struct itself (defined in
  // the spec at
  // https://arrow.apache.org/docs/format/CDataInterface.html#c.ArrowSchema.metadata)
  //
  // If metadata is non-null then it only needs to exist for the lifetime of this call,
  // a consumer who wants it to live after that must copy it to ensure lifetime.
  //
  // A producer *must not* call this concurrently from multiple different threads.
  //
  // A consumer must be prepared to receive one or more calls to this callback even
  // after calling cancel on the corresponding ArrowAsyncProducer, as cancel does not
  // guarantee it happens immediately.
  //
  // Return value: 0 if successful, `errno`-compatible error otherwise.
  //
  // If the consumer returns a non-zero return from this method, that indicates to the
  // producer that it should stop propagating data as an error occurred. After receiving
  // such a return, the only interaction with this object is for the producer to call
  // the `release` callback.
  int (*on_next_task)(struct ArrowAsyncDeviceStreamHandler* self,
                      struct ArrowAsyncTask* task, const char* metadata);

  // Handler for encountering an error. The producer should call release after
  // this returns to clean up any resources. The `code` passed in can be any error
  // code that a producer wants, but should be errno-compatible for consistency.
  //
  // If the message or metadata are non-null, they will only last as long as this
  // function call. The consumer would need to perform a copy of the data if it is
  // necessary for them to live past the lifetime of this call.
  //
  // Error metadata should be encoded as with metadata in ArrowSchema, defined in
  // the spec at
  // https://arrow.apache.org/docs/format/CDataInterface.html#c.ArrowSchema.metadata
  //
  // It is valid for this to be called by a producer with or without a preceding call
  // to ArrowAsyncProducer.request.
  //
  // This callback must not call any methods of an ArrowAsyncProducer object.
  void (*on_error)(struct ArrowAsyncDeviceStreamHandler* self, int code,
                   const char* message, const char* metadata);

  // Release callback to release any resources for the handler. Should always be
  // called by a producer when it is done utilizing a handler. No callbacks should
  // be called after this is called.
  //
  // It is valid for the release callback to be called by a producer with or without
  // a preceding call to ArrowAsyncProducer.request.
  //
  // The release callback must not call any methods of an ArrowAsyncProducer object.
  void (*release)(struct ArrowAsyncDeviceStreamHandler* self);

  // MUST be populated by the producer BEFORE calling any callbacks other than release.
  // This provides the connection between a handler and its producer, and must exist until
  // the release callback is called.
  struct ArrowAsyncProducer* producer;

  // Opaque handler-specific data
  void* private_data;
};

#endif  // ARROW_C_ASYNC_STREAM_INTERFACE

#ifdef __cplusplus
}
#endif
