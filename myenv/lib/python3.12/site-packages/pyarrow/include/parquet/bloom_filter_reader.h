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

#include "arrow/io/interfaces.h"
#include "parquet/properties.h"
#include "parquet/type_fwd.h"

namespace parquet {

class InternalFileDecryptor;
class BloomFilter;

class PARQUET_EXPORT RowGroupBloomFilterReader {
 public:
  virtual ~RowGroupBloomFilterReader() = default;

  /// \brief Read bloom filter of a column chunk.
  ///
  /// \param[in] i column ordinal of the column chunk.
  /// \returns bloom filter of the column or nullptr if it does not exist.
  /// \throws ParquetException if the index is out of bound, or read bloom
  /// filter failed.
  virtual std::unique_ptr<BloomFilter> GetColumnBloomFilter(int i) = 0;
};

/// \brief Interface for reading the bloom filter for a Parquet file.
class PARQUET_EXPORT BloomFilterReader {
 public:
  virtual ~BloomFilterReader() = default;

  /// \brief Create a BloomFilterReader instance.
  /// \returns a BloomFilterReader instance.
  /// WARNING: The returned BloomFilterReader references to all the input parameters, so
  /// it must not outlive all of the input parameters. Usually these input parameters
  /// come from the same ParquetFileReader object, so it must not outlive the reader
  /// that creates this BloomFilterReader.
  static std::unique_ptr<BloomFilterReader> Make(
      std::shared_ptr<::arrow::io::RandomAccessFile> input,
      std::shared_ptr<FileMetaData> file_metadata, const ReaderProperties& properties,
      std::shared_ptr<InternalFileDecryptor> file_decryptor = NULLPTR);

  /// \brief Get the bloom filter reader of a specific row group.
  /// \param[in] i row group ordinal to get bloom filter reader.
  /// \returns RowGroupBloomFilterReader of the specified row group. A nullptr may or may
  ///          not be returned if the bloom filter for the row group is unavailable. It
  ///          is the caller's responsibility to check the return value of follow-up calls
  ///          to the RowGroupBloomFilterReader.
  /// \throws ParquetException if the index is out of bound.
  virtual std::shared_ptr<RowGroupBloomFilterReader> RowGroup(int i) = 0;
};

}  // namespace parquet
