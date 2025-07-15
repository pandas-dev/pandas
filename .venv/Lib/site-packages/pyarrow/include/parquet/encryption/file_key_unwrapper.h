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

#include "arrow/util/concurrent_map.h"

#include "parquet/encryption/encryption.h"
#include "parquet/encryption/file_system_key_material_store.h"
#include "parquet/encryption/key_material.h"
#include "parquet/encryption/key_toolkit.h"
#include "parquet/encryption/key_toolkit_internal.h"
#include "parquet/encryption/kms_client.h"
#include "parquet/platform.h"

namespace parquet::encryption {

// This class will retrieve the key from "key metadata", following these steps:
// 1. Parse "key metadata" (see structure in KeyMetadata class).
// 2. Retrieve "key material" which can be stored inside or outside "key metadata".
// 3. Unwrap the "data encryption key" from "key material". There are 2 modes:
// 3.1. single wrapping: decrypt the wrapped "data encryption key" directly with "master
// encryption key" 3.2. double wrapping: 2 steps: 3.2.1. "key encryption key" is decrypted
// with "master encryption key" 3.2.2. "data encryption key" is decrypted with the above
// "key encryption key"
class PARQUET_EXPORT FileKeyUnwrapper : public DecryptionKeyRetriever {
 public:
  /// key_toolkit and kms_connection_config is to get KmsClient from cache or create
  /// KmsClient if it's not in the cache yet. cache_entry_lifetime_seconds is life time of
  /// KmsClient in the cache.
  /// If the file uses external key material then the Parquet file path and file
  /// system must be specified.
  FileKeyUnwrapper(std::shared_ptr<KeyToolkit> key_toolkit,
                   const KmsConnectionConfig& kms_connection_config,
                   double cache_lifetime_seconds, const std::string& file_path = "",
                   const std::shared_ptr<::arrow::fs::FileSystem>& file_system = NULLPTR);

  /// Constructor overload that takes a raw pointer to the KeyToolkit
  FileKeyUnwrapper(KeyToolkit* key_toolkit,
                   const KmsConnectionConfig& kms_connection_config,
                   double cache_lifetime_seconds, const std::string& file_path = "",
                   const std::shared_ptr<::arrow::fs::FileSystem>& file_system = NULLPTR);

  /// Constructor overload that takes a raw pointer to the KeyToolkit and
  /// accepts an existing key_material_store rather than using
  /// the file path and file system to create one when needed.
  FileKeyUnwrapper(KeyToolkit* key_toolkit,
                   const KmsConnectionConfig& kms_connection_config,
                   double cache_lifetime_seconds,
                   std::shared_ptr<FileKeyMaterialStore> key_material_store);

  /// Get the data key from key metadata
  std::string GetKey(const std::string& key_metadata) override;

  /// Get the data key along with the master key id from key material
  KeyWithMasterId GetDataEncryptionKey(const KeyMaterial& key_material);

 private:
  FileKeyUnwrapper(std::shared_ptr<KeyToolkit> key_toolkit_owner, KeyToolkit* key_toolkit,
                   const KmsConnectionConfig& kms_connection_config,
                   double cache_lifetime_seconds,
                   std::shared_ptr<FileKeyMaterialStore> key_material_store,
                   const std::string& file_path,
                   const std::shared_ptr<::arrow::fs::FileSystem>& file_system);

  std::shared_ptr<KmsClient> GetKmsClientFromConfigOrKeyMaterial(
      const KeyMaterial& key_material);

  /// A map of Key Encryption Key (KEK) ID -> KEK bytes, for the current token
  std::shared_ptr<::arrow::util::ConcurrentMap<std::string, std::string>> kek_per_kek_id_;
  std::shared_ptr<KeyToolkit> key_toolkit_owner_;
  KeyToolkit* key_toolkit_;
  KmsConnectionConfig kms_connection_config_;
  const double cache_entry_lifetime_seconds_;
  std::shared_ptr<FileKeyMaterialStore> key_material_store_;
  const std::string file_path_;
  std::shared_ptr<::arrow::fs::FileSystem> file_system_;
};

}  // namespace parquet::encryption
