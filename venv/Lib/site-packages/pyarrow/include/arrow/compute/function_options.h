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

// NOTE: API is EXPERIMENTAL and will change without going through a
// deprecation cycle.

#pragma once

#include "arrow/compute/type_fwd.h"
#include "arrow/result.h"
#include "arrow/status.h"
#include "arrow/type_fwd.h"
#include "arrow/util/visibility.h"

namespace arrow {
namespace compute {

/// \addtogroup compute-functions
/// @{

/// \brief Extension point for defining options outside libarrow (but
/// still within this project).
class ARROW_EXPORT FunctionOptionsType {
 public:
  virtual ~FunctionOptionsType() = default;

  virtual const char* type_name() const = 0;
  virtual std::string Stringify(const FunctionOptions&) const = 0;
  virtual bool Compare(const FunctionOptions&, const FunctionOptions&) const = 0;
  virtual Result<std::shared_ptr<Buffer>> Serialize(const FunctionOptions&) const;
  virtual Result<std::unique_ptr<FunctionOptions>> Deserialize(
      const Buffer& buffer) const;
  virtual std::unique_ptr<FunctionOptions> Copy(const FunctionOptions&) const = 0;
};

/// \brief Base class for specifying options configuring a function's behavior,
/// such as error handling.
class ARROW_EXPORT FunctionOptions : public util::EqualityComparable<FunctionOptions> {
 public:
  virtual ~FunctionOptions() = default;

  const FunctionOptionsType* options_type() const { return options_type_; }
  const char* type_name() const { return options_type()->type_name(); }

  bool Equals(const FunctionOptions& other) const;
  std::string ToString() const;
  std::unique_ptr<FunctionOptions> Copy() const;
  /// \brief Serialize an options struct to a buffer.
  Result<std::shared_ptr<Buffer>> Serialize() const;
  /// \brief Deserialize an options struct from a buffer.
  /// Note: this will only look for `type_name` in the default FunctionRegistry;
  /// to use a custom FunctionRegistry, look up the FunctionOptionsType, then
  /// call FunctionOptionsType::Deserialize().
  static Result<std::unique_ptr<FunctionOptions>> Deserialize(
      const std::string& type_name, const Buffer& buffer);

 protected:
  explicit FunctionOptions(const FunctionOptionsType* type) : options_type_(type) {}
  const FunctionOptionsType* options_type_;
};

ARROW_EXPORT void PrintTo(const FunctionOptions&, std::ostream*);

/// @}

}  // namespace compute
}  // namespace arrow
