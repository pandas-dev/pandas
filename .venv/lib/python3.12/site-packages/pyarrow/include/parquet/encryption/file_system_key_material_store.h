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

#include "parquet/encryption/file_key_material_store.h"
#include "parquet/exception.h"

namespace parquet::encryption {

/// A FileKeyMaterialStore that stores key material in a file system file in the same
/// folder as the Parquet file.
class PARQUET_EXPORT FileSystemKeyMaterialStore : public FileKeyMaterialStore {
 public:
  static constexpr const char kKeyMaterialFilePrefix[] = "_KEY_MATERIAL_FOR_";
  static constexpr const char kTempFilePrefix[] = "_TMP";
  static constexpr const char kKeyMaterialFileSuffix[] = ".json";

  FileSystemKeyMaterialStore() {}
  FileSystemKeyMaterialStore(std::string key_material_file_path,
                             std::shared_ptr<::arrow::fs::FileSystem> file_system);

  /// Creates a new file system key material store for a parquet file.
  /// When use_tmp_prefix is true, files are saved with an extra _TMP prefix so they don't
  /// conflict with existing external material files. This is useful during key rotation
  /// so that temporary key material files can be created while using the existing key
  /// material, before moving the key material to the non-temporary location.
  static std::shared_ptr<FileSystemKeyMaterialStore> Make(
      std::string parquet_file_path, std::shared_ptr<::arrow::fs::FileSystem> file_system,
      bool use_tmp_prefix);

  /// Add key material for one encryption key.
  void AddKeyMaterial(std::string key_id_in_file, std::string key_material) {
    key_material_map_.emplace(std::move(key_id_in_file), std::move(key_material));
  }

  /// Get key material
  std::string GetKeyMaterial(std::string key_id_in_file) {
    if (key_material_map_.empty()) {
      LoadKeyMaterialMap();
    }
    auto found = key_material_map_.find(key_id_in_file);
    if (found == key_material_map_.end()) {
      throw ParquetException("Invalid key id");
    }
    return found->second;
  }

  /// After key material was added for all keys in the given Parquet file,
  /// save material in persistent store.
  void SaveMaterial();

  /// Remove key material from persistent store. Used in key rotation.
  void RemoveMaterial();

  /// Move key material to another store. Used in key rotation.
  void MoveMaterialTo(std::shared_ptr<FileKeyMaterialStore> target_key_store);

  ///  Returns the Set of all key IDs in this store (for the given Parquet file)
  std::vector<std::string> GetKeyIDSet();

 private:
  std::string GetStorageFilePath() { return key_material_file_path_; }

  std::string BuildKeyMaterialMapJson();
  void LoadKeyMaterialMap();
  std::string key_material_file_path_;
  std::shared_ptr<::arrow::fs::FileSystem> file_system_;
  /// Maps ID of a key in Parquet file and key material
  std::unordered_map<std::string, std::string> key_material_map_;
};

}  // namespace parquet::encryption
