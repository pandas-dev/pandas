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
#include <optional>
#include <string>

#include "arrow/dataset/dataset.h"
#include "arrow/dataset/file_base.h"
#include "arrow/dataset/type_fwd.h"
#include "arrow/dataset/visibility.h"
#include "arrow/ipc/type_fwd.h"
#include "arrow/json/options.h"
#include "arrow/result.h"
#include "arrow/status.h"
#include "arrow/util/future.h"
#include "arrow/util/macros.h"

namespace arrow::dataset {

/// \addtogroup dataset-file-formats
///
/// @{

constexpr char kJsonTypeName[] = "json";

/// \brief A FileFormat implementation that reads from JSON files
class ARROW_DS_EXPORT JsonFileFormat : public FileFormat {
 public:
  JsonFileFormat();

  std::string type_name() const override { return kJsonTypeName; }

  bool Equals(const FileFormat& other) const override;

  Result<bool> IsSupported(const FileSource& source) const override;

  Result<std::shared_ptr<Schema>> Inspect(const FileSource& source) const override;

  Future<std::shared_ptr<InspectedFragment>> InspectFragment(
      const FileSource& source, const FragmentScanOptions* format_options,
      compute::ExecContext* exec_context) const override;

  Future<std::shared_ptr<FragmentScanner>> BeginScan(
      const FragmentScanRequest& scan_request, const InspectedFragment& inspected,
      const FragmentScanOptions* format_options,
      compute::ExecContext* exec_context) const override;

  Result<RecordBatchGenerator> ScanBatchesAsync(
      const std::shared_ptr<ScanOptions>& scan_options,
      const std::shared_ptr<FileFragment>& file) const override;

  Future<std::optional<int64_t>> CountRows(
      const std::shared_ptr<FileFragment>& file, compute::Expression predicate,
      const std::shared_ptr<ScanOptions>& scan_options) override;

  Result<std::shared_ptr<FileWriter>> MakeWriter(
      std::shared_ptr<io::OutputStream> destination, std::shared_ptr<Schema> schema,
      std::shared_ptr<FileWriteOptions> options,
      fs::FileLocator destination_locator) const override {
    return Status::NotImplemented("Writing JSON files is not currently supported");
  }

  std::shared_ptr<FileWriteOptions> DefaultWriteOptions() override { return NULLPTR; }
};

/// \brief Per-scan options for JSON fragments
struct ARROW_DS_EXPORT JsonFragmentScanOptions : public FragmentScanOptions {
  std::string type_name() const override { return kJsonTypeName; }

  /// @brief Options that affect JSON parsing
  ///
  /// Note: `explicit_schema` and `unexpected_field_behavior` are ignored.
  json::ParseOptions parse_options = json::ParseOptions::Defaults();

  /// @brief Options that affect JSON reading
  json::ReadOptions read_options = json::ReadOptions::Defaults();
};

/// @}

}  // namespace arrow::dataset
