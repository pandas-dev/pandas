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

#include <memory>
#include <string_view>
#include <unordered_map>

#include "arrow/result.h"
#include "arrow/util/visibility.h"

namespace arrow {
namespace json {
namespace internal {

/// This class is a helper to parse a json object from a string.
/// It uses rapidjson::Document in implementation.
class ARROW_EXPORT ObjectParser {
 public:
  ObjectParser();
  ~ObjectParser();

  Status Parse(std::string_view json);

  Result<std::string> GetString(const char* key) const;

  Result<bool> GetBool(const char* key) const;

  // Get all members of the object as a map from string keys to string values
  Result<std::unordered_map<std::string, std::string>> GetStringMap() const;

 private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace internal
}  // namespace json
}  // namespace arrow
