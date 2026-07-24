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

#include "arrow/dataset/type_fwd.h"

namespace parquet::encryption {
class CryptoFactory;
struct KmsConnectionConfig;
struct EncryptionConfiguration;
struct DecryptionConfiguration;
}  // namespace parquet::encryption

namespace arrow {
namespace dataset {

/// \brief Core configuration class encapsulating parameters for high-level encryption
/// within Parquet framework.
///
/// ParquetEncryptionConfig serves as a bridge, passing encryption-related
/// parameters to appropriate components within the Parquet library. It holds references
/// to objects defining encryption strategy, Key Management Service (KMS) configuration,
/// and specific encryption configurations for Parquet data.
struct ARROW_DS_EXPORT ParquetEncryptionConfig {
  ///  Shared pointer to CryptoFactory object, responsible for creating cryptographic
  ///  components like encryptors and decryptors.
  std::shared_ptr<parquet::encryption::CryptoFactory> crypto_factory;

  ///  Shared pointer to KmsConnectionConfig object, holding configuration parameters for
  ///  connecting to a Key Management Service (KMS).
  std::shared_ptr<parquet::encryption::KmsConnectionConfig> kms_connection_config;

  ///  Shared pointer to EncryptionConfiguration object, defining specific encryption
  ///  settings for Parquet data, like keys for different columns.
  std::shared_ptr<parquet::encryption::EncryptionConfiguration> encryption_config;
};

/// \brief Core configuration class encapsulating parameters for high-level decryption
/// within Parquet framework.
///
/// ParquetDecryptionConfig is designed to pass decryption-related parameters to
/// appropriate decryption components within Parquet library. It holds references to
/// objects defining decryption strategy, Key Management Service (KMS) configuration,
/// and specific decryption configurations for reading encrypted Parquet data.
struct ARROW_DS_EXPORT ParquetDecryptionConfig {
  ///  Shared pointer to CryptoFactory object, pivotal in creating cryptographic
  ///  components for decryption process.
  std::shared_ptr<parquet::encryption::CryptoFactory> crypto_factory;

  ///  Shared pointer to KmsConnectionConfig object, containing parameters for connecting
  ///  to a Key Management Service (KMS) during decryption.
  std::shared_ptr<parquet::encryption::KmsConnectionConfig> kms_connection_config;

  ///  Shared pointer to DecryptionConfiguration object, specifying decryption settings
  ///  for reading encrypted Parquet data.
  std::shared_ptr<parquet::encryption::DecryptionConfiguration> decryption_config;
};

}  // namespace dataset
}  // namespace arrow
