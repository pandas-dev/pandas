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

#include <bitset>
#include <string>

#include "arrow/type_fwd.h"
#include "arrow/util/visibility.h"

namespace arrow {

ARROW_EXPORT
const char* DeviceAllocationTypeToCStr(DeviceAllocationType type);

class ARROW_EXPORT DeviceAllocationTypeSet {
 private:
  std::bitset<kDeviceAllocationTypeMax + 1> device_type_bitset_;

 public:
  /// \brief Construct an empty set of device types.
  DeviceAllocationTypeSet() = default;

  /// \brief Construct a set of device types with a single device type.
  DeviceAllocationTypeSet(  // NOLINT implicit construction
      DeviceAllocationType accepted_device_type) {
    add(accepted_device_type);
  }

  /// \brief Construct a set of device types containing only "kCPU".
  static DeviceAllocationTypeSet CpuOnly() {
    return DeviceAllocationTypeSet{DeviceAllocationType::kCPU};
  }

  /// \brief Construct a set of device types containing all device types.
  static DeviceAllocationTypeSet All() {
    DeviceAllocationTypeSet all;
    all.device_type_bitset_.set();
    // Don't set the invalid enum values.
    all.device_type_bitset_.reset(0);
    all.device_type_bitset_.reset(5);
    all.device_type_bitset_.reset(6);
    return all;
  }

  /// \brief Add a device type to the set of device types.
  void add(DeviceAllocationType device_type) {
    device_type_bitset_.set(static_cast<int>(device_type));
  }

  /// \brief Remove a device type from the set of device types.
  void remove(DeviceAllocationType device_type) {
    device_type_bitset_.reset(static_cast<int>(device_type));
  }

  /// \brief Return true iff the set only contains the CPU device type.
  bool is_cpu_only() const {
    return device_type_bitset_ == CpuOnly().device_type_bitset_;
  }

  /// \brief Return true if the set of accepted device types includes the
  /// device type.
  bool contains(DeviceAllocationType device_type) const {
    return device_type_bitset_.test(static_cast<int>(device_type));
  }

  /// \brief Add all device types from another set to this set.
  void Add(DeviceAllocationTypeSet other) {
    device_type_bitset_ |= other.device_type_bitset_;
  }

  /// \brief Return true if the set of accepted device types includes all the
  /// device types in the other set.
  bool Contains(DeviceAllocationTypeSet other) const {
    // other \subseteq this <==> (other \intersect this == other)
    return (other.device_type_bitset_ & device_type_bitset_) == other.device_type_bitset_;
  }

  std::string ToString() const;
};

}  // namespace arrow
