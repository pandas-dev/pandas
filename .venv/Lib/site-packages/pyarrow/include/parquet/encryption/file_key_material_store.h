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

#include <set>
#include <string>
#include <unordered_map>

#include "arrow/filesystem/filesystem.h"
#include "parquet/platform.h"

namespace parquet::encryption {

/// Stores encryption key material outside the Parquet file, for example in a separate
/// small file in the same folder. This is important for “key rotation”, when MEKs have to
/// be changed (if compromised; or periodically, just in case) - without modifying the
/// Parquet files (often  immutable).
class PARQUET_EXPORT FileKeyMaterialStore {
 public:
  /// Add key material for one encryption key.
  virtual void AddKeyMaterial(std::string key_id_in_file, std::string key_material) = 0;

  /// Get key material
  virtual std::string GetKeyMaterial(std::string key_id_in_file) = 0;

  /// After key material was added for all keys in the given Parquet file,
  /// save material in persistent store.
  virtual void SaveMaterial() = 0;

  /// Remove key material from persistent store. Used in key rotation.
  virtual void RemoveMaterial() = 0;

  /// Move key material to another store. Used in key rotation.
  virtual void MoveMaterialTo(std::shared_ptr<FileKeyMaterialStore> target_key_store) = 0;

  /// Returns the Set of all key IDs in this store (for the given Parquet file)
  virtual std::vector<std::string> GetKeyIDSet() = 0;

  virtual ~FileKeyMaterialStore() {}
};

}  // namespace parquet::encryption
