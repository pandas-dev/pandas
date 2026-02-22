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

#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>

#include "arrow/result.h"
#include "arrow/status.h"
#include "arrow/type_fwd.h"
#include "arrow/util/macros.h"
#include "arrow/util/visibility.h"

namespace arrow {

namespace internal {

///////////////////////////////////////////////////////////////////////
// Helper tracking memory statistics

/// \brief Memory pool statistics
///
/// 64-byte aligned so that all atomic values are on the same cache line.
class alignas(64) MemoryPoolStats {
 private:
  // All atomics are updated according to Acquire-Release ordering.
  // https://en.cppreference.com/w/cpp/atomic/memory_order#Release-Acquire_ordering
  //
  // max_memory_, total_allocated_bytes_, and num_allocs_ only go up (they are
  // monotonically increasing) which can allow some optimizations.
  std::atomic<int64_t> max_memory_{0};
  std::atomic<int64_t> bytes_allocated_{0};
  std::atomic<int64_t> total_allocated_bytes_{0};
  std::atomic<int64_t> num_allocs_{0};

 public:
  int64_t max_memory() const { return max_memory_.load(std::memory_order_acquire); }

  int64_t bytes_allocated() const {
    return bytes_allocated_.load(std::memory_order_acquire);
  }

  int64_t total_bytes_allocated() const {
    return total_allocated_bytes_.load(std::memory_order_acquire);
  }

  int64_t num_allocations() const { return num_allocs_.load(std::memory_order_acquire); }

  inline void DidAllocateBytes(int64_t size) {
    // Issue the load before everything else. max_memory_ is monotonically increasing,
    // so we can use a relaxed load before the read-modify-write.
    auto max_memory = max_memory_.load(std::memory_order_relaxed);
    const auto old_bytes_allocated =
        bytes_allocated_.fetch_add(size, std::memory_order_acq_rel);
    // Issue store operations on values that we don't depend on to proceed
    // with execution. When done, max_memory and old_bytes_allocated have
    // a higher chance of being available on CPU registers. This also has the
    // nice side-effect of putting 3 atomic stores close to each other in the
    // instruction stream.
    total_allocated_bytes_.fetch_add(size, std::memory_order_acq_rel);
    num_allocs_.fetch_add(1, std::memory_order_acq_rel);

    // If other threads are updating max_memory_ concurrently we leave the loop without
    // updating knowing that it already reached a value even higher than ours.
    const auto allocated = old_bytes_allocated + size;
    while (max_memory < allocated && !max_memory_.compare_exchange_weak(
                                         /*expected=*/max_memory, /*desired=*/allocated,
                                         std::memory_order_acq_rel)) {
    }
  }

  inline void DidReallocateBytes(int64_t old_size, int64_t new_size) {
    if (new_size > old_size) {
      DidAllocateBytes(new_size - old_size);
    } else {
      DidFreeBytes(old_size - new_size);
    }
  }

  inline void DidFreeBytes(int64_t size) {
    bytes_allocated_.fetch_sub(size, std::memory_order_acq_rel);
  }
};

}  // namespace internal

/// Base class for memory allocation on the CPU.
///
/// Besides tracking the number of allocated bytes, the allocator also should
/// take care of the required 64-byte alignment.
class ARROW_EXPORT MemoryPool {
 public:
  virtual ~MemoryPool() = default;

  /// \brief EXPERIMENTAL. Create a new instance of the default MemoryPool
  static std::unique_ptr<MemoryPool> CreateDefault();

  /// Allocate a new memory region of at least size bytes.
  ///
  /// The allocated region shall be 64-byte aligned.
  Status Allocate(int64_t size, uint8_t** out) {
    return Allocate(size, kDefaultBufferAlignment, out);
  }

  /// Allocate a new memory region of at least size bytes aligned to alignment.
  virtual Status Allocate(int64_t size, int64_t alignment, uint8_t** out) = 0;

  /// Resize an already allocated memory section.
  ///
  /// As by default most default allocators on a platform don't support aligned
  /// reallocation, this function can involve a copy of the underlying data.
  virtual Status Reallocate(int64_t old_size, int64_t new_size, int64_t alignment,
                            uint8_t** ptr) = 0;
  Status Reallocate(int64_t old_size, int64_t new_size, uint8_t** ptr) {
    return Reallocate(old_size, new_size, kDefaultBufferAlignment, ptr);
  }

  /// Free an allocated region.
  ///
  /// @param buffer Pointer to the start of the allocated memory region
  /// @param size Allocated size located at buffer. An allocator implementation
  ///   may use this for tracking the amount of allocated bytes as well as for
  ///   faster deallocation if supported by its backend.
  /// @param alignment The alignment of the allocation. Defaults to 64 bytes.
  virtual void Free(uint8_t* buffer, int64_t size, int64_t alignment) = 0;
  void Free(uint8_t* buffer, int64_t size) {
    Free(buffer, size, kDefaultBufferAlignment);
  }

  /// Return unused memory to the OS
  ///
  /// Only applies to allocators that hold onto unused memory.  This will be
  /// best effort, a memory pool may not implement this feature or may be
  /// unable to fulfill the request due to fragmentation.
  virtual void ReleaseUnused() {}

  /// Print statistics
  ///
  /// Print allocation statistics on stderr. The output format is
  /// implementation-specific. Not all memory pools implement this method.
  virtual void PrintStats() {}

  /// The number of bytes that were allocated and not yet free'd through
  /// this allocator.
  virtual int64_t bytes_allocated() const = 0;

  /// Return peak memory allocation in this memory pool
  ///
  /// \return Maximum bytes allocated. If not known (or not implemented),
  /// returns -1
  virtual int64_t max_memory() const;

  /// The number of bytes that were allocated.
  virtual int64_t total_bytes_allocated() const = 0;

  /// The number of allocations or reallocations that were requested.
  virtual int64_t num_allocations() const = 0;

  /// The name of the backend used by this MemoryPool (e.g. "system" or "jemalloc").
  virtual std::string backend_name() const = 0;

 protected:
  MemoryPool() = default;
};

class ARROW_EXPORT LoggingMemoryPool : public MemoryPool {
 public:
  explicit LoggingMemoryPool(MemoryPool* pool);
  ~LoggingMemoryPool() override = default;

  using MemoryPool::Allocate;
  using MemoryPool::Free;
  using MemoryPool::Reallocate;

  Status Allocate(int64_t size, int64_t alignment, uint8_t** out) override;
  Status Reallocate(int64_t old_size, int64_t new_size, int64_t alignment,
                    uint8_t** ptr) override;
  void Free(uint8_t* buffer, int64_t size, int64_t alignment) override;
  void ReleaseUnused() override;
  void PrintStats() override;

  int64_t bytes_allocated() const override;

  int64_t max_memory() const override;

  int64_t total_bytes_allocated() const override;

  int64_t num_allocations() const override;

  std::string backend_name() const override;

 private:
  MemoryPool* pool_;
};

/// Derived class for memory allocation.
///
/// Tracks the number of bytes and maximum memory allocated through its direct
/// calls. Actual allocation is delegated to MemoryPool class.
class ARROW_EXPORT ProxyMemoryPool : public MemoryPool {
 public:
  explicit ProxyMemoryPool(MemoryPool* pool);
  ~ProxyMemoryPool() override;

  using MemoryPool::Allocate;
  using MemoryPool::Free;
  using MemoryPool::Reallocate;

  Status Allocate(int64_t size, int64_t alignment, uint8_t** out) override;
  Status Reallocate(int64_t old_size, int64_t new_size, int64_t alignment,
                    uint8_t** ptr) override;
  void Free(uint8_t* buffer, int64_t size, int64_t alignment) override;
  void ReleaseUnused() override;
  void PrintStats() override;

  int64_t bytes_allocated() const override;

  int64_t max_memory() const override;

  int64_t total_bytes_allocated() const override;

  int64_t num_allocations() const override;

  std::string backend_name() const override;

 private:
  class ProxyMemoryPoolImpl;
  std::unique_ptr<ProxyMemoryPoolImpl> impl_;
};

/// EXPERIMENTAL MemoryPool wrapper with an upper limit
///
/// Checking for limits is not done in a fully thread-safe way, therefore
/// multi-threaded allocations might be able to go successfully above the
/// configured limit.
class ARROW_EXPORT CappedMemoryPool : public MemoryPool {
 public:
  CappedMemoryPool(MemoryPool* wrapped_pool, int64_t bytes_allocated_limit)
      : wrapped_(wrapped_pool), bytes_allocated_limit_(bytes_allocated_limit) {}

  using MemoryPool::Allocate;
  using MemoryPool::Reallocate;

  Status Allocate(int64_t size, int64_t alignment, uint8_t** out) override;
  Status Reallocate(int64_t old_size, int64_t new_size, int64_t alignment,
                    uint8_t** ptr) override;
  void Free(uint8_t* buffer, int64_t size, int64_t alignment) override;

  void ReleaseUnused() override { wrapped_->ReleaseUnused(); }

  void PrintStats() override { wrapped_->PrintStats(); }

  int64_t bytes_allocated() const override { return wrapped_->bytes_allocated(); }

  int64_t max_memory() const override { return wrapped_->max_memory(); }

  int64_t total_bytes_allocated() const override {
    return wrapped_->total_bytes_allocated();
  }

  int64_t num_allocations() const override { return wrapped_->num_allocations(); }

  std::string backend_name() const override { return wrapped_->backend_name(); }

 private:
  Status OutOfMemory(int64_t current_allocated, int64_t requested) const;

  MemoryPool* wrapped_;
  const int64_t bytes_allocated_limit_;
};

/// \brief Return a process-wide memory pool based on the system allocator.
ARROW_EXPORT MemoryPool* system_memory_pool();

/// \brief Return a process-wide memory pool based on jemalloc.
///
/// May return NotImplemented if jemalloc is not available.
ARROW_EXPORT Status jemalloc_memory_pool(MemoryPool** out);

/// \brief Set jemalloc memory page purging behavior for future-created arenas
/// to the indicated number of milliseconds. See dirty_decay_ms and
/// muzzy_decay_ms options in jemalloc for a description of what these do. The
/// default is configured to 1000 (1 second) which releases memory more
/// aggressively to the operating system than the jemalloc default of 10
/// seconds. If you set the value to 0, dirty / muzzy pages will be released
/// immediately rather than with a time decay, but this may reduce application
/// performance.
ARROW_EXPORT
Status jemalloc_set_decay_ms(int ms);

/// \brief Get basic statistics from jemalloc's mallctl.
/// See the MALLCTL NAMESPACE section in jemalloc project documentation for
/// available stats.
ARROW_EXPORT
Result<int64_t> jemalloc_get_stat(const char* name);

/// \brief Reset the counter for peak bytes allocated in the calling thread to zero.
/// This affects subsequent calls to thread.peak.read, but not the values returned by
/// thread.allocated or thread.deallocated.
ARROW_EXPORT
Status jemalloc_peak_reset();

/// \brief Print summary statistics in human-readable form to stderr.
/// See malloc_stats_print documentation in jemalloc project documentation for
/// available opt flags.
ARROW_EXPORT
Status jemalloc_stats_print(const char* opts = "");

/// \brief Print summary statistics in human-readable form using a callback
/// See malloc_stats_print documentation in jemalloc project documentation for
/// available opt flags.
ARROW_EXPORT
Status jemalloc_stats_print(std::function<void(const char*)> write_cb,
                            const char* opts = "");

/// \brief Get summary statistics in human-readable form.
/// See malloc_stats_print documentation in jemalloc project documentation for
/// available opt flags.
ARROW_EXPORT
Result<std::string> jemalloc_stats_string(const char* opts = "");

/// \brief Return a process-wide memory pool based on mimalloc.
///
/// May return NotImplemented if mimalloc is not available.
ARROW_EXPORT Status mimalloc_memory_pool(MemoryPool** out);

/// \brief Return the names of the backends supported by this Arrow build.
ARROW_EXPORT std::vector<std::string> SupportedMemoryBackendNames();

}  // namespace arrow
