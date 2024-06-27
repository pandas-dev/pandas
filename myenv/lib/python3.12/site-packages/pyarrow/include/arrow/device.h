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

#include <cstdint>
#include <functional>
#include <memory>
#include <string>

#include "arrow/io/type_fwd.h"
#include "arrow/result.h"
#include "arrow/status.h"
#include "arrow/type_fwd.h"
#include "arrow/util/compare.h"
#include "arrow/util/macros.h"
#include "arrow/util/visibility.h"

namespace arrow {

/// \brief EXPERIMENTAL: Device type enum which matches up with C Data Device types
enum class DeviceAllocationType : char {
  kCPU = 1,
  kCUDA = 2,
  kCUDA_HOST = 3,
  kOPENCL = 4,
  kVULKAN = 7,
  kMETAL = 8,
  kVPI = 9,
  kROCM = 10,
  kROCM_HOST = 11,
  kEXT_DEV = 12,
  kCUDA_MANAGED = 13,
  kONEAPI = 14,
  kWEBGPU = 15,
  kHEXAGON = 16,
};

class MemoryManager;

/// \brief EXPERIMENTAL: Abstract interface for hardware devices
///
/// This object represents a device with access to some memory spaces.
/// When handling a Buffer or raw memory address, it allows deciding in which
/// context the raw memory address should be interpreted
/// (e.g. CPU-accessible memory, or embedded memory on some particular GPU).
class ARROW_EXPORT Device : public std::enable_shared_from_this<Device>,
                            public util::EqualityComparable<Device> {
 public:
  virtual ~Device();

  /// \brief A shorthand for this device's type.
  ///
  /// The returned value is different for each device class, but is the
  /// same for all instances of a given class.  It can be used as a replacement
  /// for RTTI.
  virtual const char* type_name() const = 0;

  /// \brief A human-readable description of the device.
  ///
  /// The returned value should be detailed enough to distinguish between
  /// different instances, where necessary.
  virtual std::string ToString() const = 0;

  /// \brief Whether this instance points to the same device as another one.
  virtual bool Equals(const Device&) const = 0;

  /// \brief A device ID to identify this device if there are multiple of this type.
  ///
  /// If there is no "device_id" equivalent (such as for the main CPU device on
  /// non-numa systems) returns -1.
  virtual int64_t device_id() const { return -1; }

  /// \brief Whether this device is the main CPU device.
  ///
  /// This shorthand method is very useful when deciding whether a memory address
  /// is CPU-accessible.
  bool is_cpu() const { return is_cpu_; }

  /// \brief Return a MemoryManager instance tied to this device
  ///
  /// The returned instance uses default parameters for this device type's
  /// MemoryManager implementation.  Some devices also allow constructing
  /// MemoryManager instances with non-default parameters.
  virtual std::shared_ptr<MemoryManager> default_memory_manager() = 0;

  /// \brief Return the DeviceAllocationType of this device
  virtual DeviceAllocationType device_type() const = 0;

  class SyncEvent;

  /// \brief EXPERIMENTAL: An opaque wrapper for Device-specific streams
  ///
  /// In essence this is just a wrapper around a void* to represent the
  /// standard concept of a stream/queue on a device. Derived classes
  /// should be trivially constructible from it's device-specific counterparts.
  class ARROW_EXPORT Stream {
   public:
    using release_fn_t = std::function<void(void*)>;

    virtual ~Stream() = default;

    virtual const void* get_raw() const { return stream_.get(); }

    /// \brief Make the stream wait on the provided event.
    ///
    /// Tells the stream that it should wait until the synchronization
    /// event is completed without blocking the CPU.
    virtual Status WaitEvent(const SyncEvent&) = 0;

    /// \brief Blocks the current thread until a stream's remaining tasks are completed
    virtual Status Synchronize() const = 0;

   protected:
    explicit Stream(void* stream, release_fn_t release_stream)
        : stream_{stream, release_stream} {}

    std::unique_ptr<void, release_fn_t> stream_;
  };

  virtual Result<std::shared_ptr<Stream>> MakeStream() { return NULLPTR; }

  /// \brief Create a new device stream
  ///
  /// This should create the appropriate stream type for the device,
  /// derived from Device::Stream to allow for stream ordered events
  /// and memory allocations.
  virtual Result<std::shared_ptr<Stream>> MakeStream(unsigned int flags) {
    return NULLPTR;
  }

  /// @brief Wrap an existing device stream alongside a release function
  ///
  /// @param device_stream a pointer to the stream to wrap
  /// @param release_fn a function to call during destruction, `nullptr` or
  ///        a no-op function can be passed to indicate ownership is maintained
  ///        externally
  virtual Result<std::shared_ptr<Stream>> WrapStream(void* device_stream,
                                                     Stream::release_fn_t release_fn) {
    return NULLPTR;
  }

  /// \brief EXPERIMENTAL: An object that provides event/stream sync primitives
  class ARROW_EXPORT SyncEvent {
   public:
    using release_fn_t = std::function<void(void*)>;

    virtual ~SyncEvent() = default;

    void* get_raw() { return sync_event_.get(); }

    /// @brief Block until sync event is completed.
    virtual Status Wait() = 0;

    /// @brief Record the wrapped event on the stream so it triggers
    /// the event when the stream gets to that point in its queue.
    virtual Status Record(const Stream&) = 0;

   protected:
    /// If creating this with a passed in event, the caller must ensure
    /// that the event lives until clear_event is called on this as it
    /// won't own it.
    explicit SyncEvent(void* sync_event, release_fn_t release_sync_event)
        : sync_event_{sync_event, release_sync_event} {}

    std::unique_ptr<void, release_fn_t> sync_event_;
  };

 protected:
  ARROW_DISALLOW_COPY_AND_ASSIGN(Device);
  explicit Device(bool is_cpu = false) : is_cpu_(is_cpu) {}

  bool is_cpu_;
};

/// \brief EXPERIMENTAL: An object that provides memory management primitives
///
/// A MemoryManager is always tied to a particular Device instance.
/// It can also have additional parameters (such as a MemoryPool to
/// allocate CPU memory).
class ARROW_EXPORT MemoryManager : public std::enable_shared_from_this<MemoryManager> {
 public:
  virtual ~MemoryManager();

  /// \brief The device this MemoryManager is tied to
  const std::shared_ptr<Device>& device() const { return device_; }

  /// \brief Whether this MemoryManager is tied to the main CPU device.
  ///
  /// This shorthand method is very useful when deciding whether a memory address
  /// is CPU-accessible.
  bool is_cpu() const { return device_->is_cpu(); }

  /// \brief Create a RandomAccessFile to read a particular buffer.
  ///
  /// The given buffer must be tied to this MemoryManager.
  ///
  /// See also the Buffer::GetReader shorthand.
  virtual Result<std::shared_ptr<io::RandomAccessFile>> GetBufferReader(
      std::shared_ptr<Buffer> buf) = 0;

  /// \brief Create a OutputStream to write to a particular buffer.
  ///
  /// The given buffer must be mutable and tied to this MemoryManager.
  /// The returned stream object writes into the buffer's underlying memory
  /// (but it won't resize it).
  ///
  /// See also the Buffer::GetWriter shorthand.
  virtual Result<std::shared_ptr<io::OutputStream>> GetBufferWriter(
      std::shared_ptr<Buffer> buf) = 0;

  /// \brief Allocate a (mutable) Buffer
  ///
  /// The buffer will be allocated in the device's memory.
  virtual Result<std::unique_ptr<Buffer>> AllocateBuffer(int64_t size) = 0;

  /// \brief Copy a Buffer to a destination MemoryManager
  ///
  /// See also the Buffer::Copy shorthand.
  static Result<std::shared_ptr<Buffer>> CopyBuffer(
      const std::shared_ptr<Buffer>& source, const std::shared_ptr<MemoryManager>& to);

  /// \brief Copy a non-owned Buffer to a destination MemoryManager
  ///
  /// This is useful for cases where the source memory area is externally managed
  /// (its lifetime not tied to the source Buffer), otherwise please use CopyBuffer().
  static Result<std::unique_ptr<Buffer>> CopyNonOwned(
      const Buffer& source, const std::shared_ptr<MemoryManager>& to);

  /// \brief Make a no-copy Buffer view in a destination MemoryManager
  ///
  /// See also the Buffer::View shorthand.
  static Result<std::shared_ptr<Buffer>> ViewBuffer(
      const std::shared_ptr<Buffer>& source, const std::shared_ptr<MemoryManager>& to);

  /// \brief Create a new SyncEvent.
  ///
  /// This version should construct the appropriate event for the device and
  /// provide the unique_ptr with the correct deleter for the event type.
  /// If the device does not require or work with any synchronization, it is
  /// allowed for it to return a nullptr.
  virtual Result<std::shared_ptr<Device::SyncEvent>> MakeDeviceSyncEvent();

  /// \brief Wrap an event into a SyncEvent.
  ///
  /// @param sync_event passed in sync_event (should be a pointer to the appropriate type)
  /// @param release_sync_event destructor to free sync_event. `nullptr` may be
  ///        passed to indicate that no destruction/freeing is necessary
  virtual Result<std::shared_ptr<Device::SyncEvent>> WrapDeviceSyncEvent(
      void* sync_event, Device::SyncEvent::release_fn_t release_sync_event);

 protected:
  ARROW_DISALLOW_COPY_AND_ASSIGN(MemoryManager);

  explicit MemoryManager(const std::shared_ptr<Device>& device) : device_(device) {}

  // Default implementations always return nullptr, should be overridden
  // by subclasses that support data transfer.
  // (returning nullptr means unsupported copy / view)
  // In CopyBufferFrom and ViewBufferFrom, the `from` parameter is guaranteed to
  // be equal to `buf->memory_manager()`.
  virtual Result<std::shared_ptr<Buffer>> CopyBufferFrom(
      const std::shared_ptr<Buffer>& buf, const std::shared_ptr<MemoryManager>& from);
  virtual Result<std::shared_ptr<Buffer>> CopyBufferTo(
      const std::shared_ptr<Buffer>& buf, const std::shared_ptr<MemoryManager>& to);
  virtual Result<std::unique_ptr<Buffer>> CopyNonOwnedFrom(
      const Buffer& buf, const std::shared_ptr<MemoryManager>& from);
  virtual Result<std::unique_ptr<Buffer>> CopyNonOwnedTo(
      const Buffer& buf, const std::shared_ptr<MemoryManager>& to);
  virtual Result<std::shared_ptr<Buffer>> ViewBufferFrom(
      const std::shared_ptr<Buffer>& buf, const std::shared_ptr<MemoryManager>& from);
  virtual Result<std::shared_ptr<Buffer>> ViewBufferTo(
      const std::shared_ptr<Buffer>& buf, const std::shared_ptr<MemoryManager>& to);

  std::shared_ptr<Device> device_;
};

// ----------------------------------------------------------------------
// CPU backend implementation

class ARROW_EXPORT CPUDevice : public Device {
 public:
  const char* type_name() const override;
  std::string ToString() const override;
  bool Equals(const Device&) const override;
  DeviceAllocationType device_type() const override { return DeviceAllocationType::kCPU; }

  std::shared_ptr<MemoryManager> default_memory_manager() override;

  /// \brief Return the global CPUDevice instance
  static std::shared_ptr<Device> Instance();

  /// \brief Create a MemoryManager
  ///
  /// The returned MemoryManager will use the given MemoryPool for allocations.
  static std::shared_ptr<MemoryManager> memory_manager(MemoryPool* pool);

 protected:
  CPUDevice() : Device(true) {}
};

class ARROW_EXPORT CPUMemoryManager : public MemoryManager {
 public:
  Result<std::shared_ptr<io::RandomAccessFile>> GetBufferReader(
      std::shared_ptr<Buffer> buf) override;
  Result<std::shared_ptr<io::OutputStream>> GetBufferWriter(
      std::shared_ptr<Buffer> buf) override;

  Result<std::unique_ptr<Buffer>> AllocateBuffer(int64_t size) override;

  /// \brief Return the MemoryPool associated with this MemoryManager.
  MemoryPool* pool() const { return pool_; }

 protected:
  CPUMemoryManager(const std::shared_ptr<Device>& device, MemoryPool* pool)
      : MemoryManager(device), pool_(pool) {}

  static std::shared_ptr<MemoryManager> Make(const std::shared_ptr<Device>& device,
                                             MemoryPool* pool = default_memory_pool());

  Result<std::shared_ptr<Buffer>> CopyBufferFrom(
      const std::shared_ptr<Buffer>& buf,
      const std::shared_ptr<MemoryManager>& from) override;
  Result<std::shared_ptr<Buffer>> CopyBufferTo(
      const std::shared_ptr<Buffer>& buf,
      const std::shared_ptr<MemoryManager>& to) override;
  Result<std::unique_ptr<Buffer>> CopyNonOwnedFrom(
      const Buffer& buf, const std::shared_ptr<MemoryManager>& from) override;
  Result<std::unique_ptr<Buffer>> CopyNonOwnedTo(
      const Buffer& buf, const std::shared_ptr<MemoryManager>& to) override;
  Result<std::shared_ptr<Buffer>> ViewBufferFrom(
      const std::shared_ptr<Buffer>& buf,
      const std::shared_ptr<MemoryManager>& from) override;
  Result<std::shared_ptr<Buffer>> ViewBufferTo(
      const std::shared_ptr<Buffer>& buf,
      const std::shared_ptr<MemoryManager>& to) override;

  MemoryPool* pool_;

  friend std::shared_ptr<MemoryManager> CPUDevice::memory_manager(MemoryPool* pool);
  ARROW_FRIEND_EXPORT friend std::shared_ptr<MemoryManager> default_cpu_memory_manager();
};

/// \brief Return the default CPU MemoryManager instance
///
/// The returned singleton instance uses the default MemoryPool.
/// This function is a faster spelling of
/// `CPUDevice::Instance()->default_memory_manager()`.
ARROW_EXPORT
std::shared_ptr<MemoryManager> default_cpu_memory_manager();

using DeviceMapper =
    std::function<Result<std::shared_ptr<MemoryManager>>(int64_t device_id)>;

/// \brief Register a function to retrieve a MemoryManager for a Device type
///
/// This registers the device type globally. A specific device type can only
/// be registered once. This method is thread-safe.
///
/// Currently, this registry is only used for importing data through the C Device
/// Data Interface (for the default Device to MemoryManager mapper in
/// arrow::ImportDeviceArray/ImportDeviceRecordBatch).
///
/// \param[in] device_type the device type for which to register a MemoryManager
/// \param[in] mapper function that takes a device id and returns the appropriate
/// MemoryManager for the registered device type and given device id
/// \return Status
ARROW_EXPORT
Status RegisterDeviceMapper(DeviceAllocationType device_type, DeviceMapper mapper);

/// \brief Get the registered function to retrieve a MemoryManager for the
/// given Device type
///
/// \param[in] device_type the device type
/// \return function that takes a device id and returns the appropriate
/// MemoryManager for the registered device type and given device id
ARROW_EXPORT
Result<DeviceMapper> GetDeviceMapper(DeviceAllocationType device_type);

}  // namespace arrow
