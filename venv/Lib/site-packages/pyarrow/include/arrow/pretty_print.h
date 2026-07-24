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

#include <iosfwd>
#include <string>
#include <utility>

#include "arrow/util/visibility.h"

namespace arrow {

class Array;
class ChunkedArray;
class RecordBatch;
class Schema;
class Status;
class Table;

/// \class PrettyPrintDelimiters
/// \brief Options for controlling which delimiters to use when printing
/// an Array or ChunkedArray.
struct ARROW_EXPORT PrettyPrintDelimiters {
  /// Delimiter to use when opening an Array or ChunkedArray (e.g. "[")
  std::string open = "[";

  /// Delimiter to use when closing an Array or ChunkedArray (e.g. "]")
  std::string close = "]";

  /// Delimiter for separating individual elements of an Array (e.g. ","),
  /// or individual chunks of a ChunkedArray
  std::string element = ",";

  /// Create a PrettyPrintDelimiters instance with default values
  static PrettyPrintDelimiters Defaults() { return PrettyPrintDelimiters(); }
};

/// \class PrettyPrintOptions
/// \brief Options for controlling how various Arrow types should be printed.
struct ARROW_EXPORT PrettyPrintOptions {
  PrettyPrintOptions() = default;

  PrettyPrintOptions(int indent,  // NOLINT runtime/explicit
                     int window = 10, int indent_size = 2, std::string null_rep = "null",
                     bool skip_new_lines = false, bool truncate_metadata = true,
                     int container_window = 2, int element_size_limit = 100)
      : indent(indent),
        indent_size(indent_size),
        window(window),
        container_window(container_window),
        null_rep(std::move(null_rep)),
        skip_new_lines(skip_new_lines),
        truncate_metadata(truncate_metadata),
        element_size_limit(element_size_limit) {}

  /// Create a PrettyPrintOptions instance with default values
  static PrettyPrintOptions Defaults() { return PrettyPrintOptions(); }

  /// Number of spaces to shift entire formatted object to the right
  int indent = 0;

  /// Size of internal indents
  int indent_size = 2;

  /// Maximum number of elements to show at the beginning and at the end.
  int window = 10;

  /// Maximum number of elements to show at the beginning and at the end, for elements
  /// that are containers (that is, list in ListArray and chunks in ChunkedArray)
  int container_window = 2;

  /// String to use for representing a null value, defaults to "null"
  std::string null_rep = "null";

  /// Skip new lines between elements, defaults to false
  bool skip_new_lines = false;

  /// Limit display of each KeyValueMetadata key/value pair to a single line at
  /// 80 character width
  bool truncate_metadata = true;

  /// If true, display field metadata when pretty-printing a Schema
  bool show_field_metadata = true;

  /// If true, display schema metadata when pretty-printing a Schema
  bool show_schema_metadata = true;

  /// Limit each element to specified number of characters, defaults to 100
  int element_size_limit = 100;

  /// Delimiters to use when printing an Array
  PrettyPrintDelimiters array_delimiters = PrettyPrintDelimiters::Defaults();

  /// Delimiters to use when printing a ChunkedArray
  PrettyPrintDelimiters chunked_array_delimiters = PrettyPrintDelimiters::Defaults();
};

/// \brief Print human-readable representation of RecordBatch
ARROW_EXPORT
Status PrettyPrint(const RecordBatch& batch, int indent, std::ostream* sink);

ARROW_EXPORT
Status PrettyPrint(const RecordBatch& batch, const PrettyPrintOptions& options,
                   std::ostream* sink);

/// \brief Print human-readable representation of Table
ARROW_EXPORT
Status PrettyPrint(const Table& table, const PrettyPrintOptions& options,
                   std::ostream* sink);

/// \brief Print human-readable representation of Array
ARROW_EXPORT
Status PrettyPrint(const Array& arr, int indent, std::ostream* sink);

/// \brief Print human-readable representation of Array
ARROW_EXPORT
Status PrettyPrint(const Array& arr, const PrettyPrintOptions& options,
                   std::ostream* sink);

/// \brief Print human-readable representation of Array
ARROW_EXPORT
Status PrettyPrint(const Array& arr, const PrettyPrintOptions& options,
                   std::string* result);

/// \brief Print human-readable representation of ChunkedArray
ARROW_EXPORT
Status PrettyPrint(const ChunkedArray& chunked_arr, const PrettyPrintOptions& options,
                   std::ostream* sink);

/// \brief Print human-readable representation of ChunkedArray
ARROW_EXPORT
Status PrettyPrint(const ChunkedArray& chunked_arr, const PrettyPrintOptions& options,
                   std::string* result);

ARROW_EXPORT
Status PrettyPrint(const Schema& schema, const PrettyPrintOptions& options,
                   std::ostream* sink);

ARROW_EXPORT
Status PrettyPrint(const Schema& schema, const PrettyPrintOptions& options,
                   std::string* result);

ARROW_EXPORT
Status DebugPrint(const Array& arr, int indent);

}  // namespace arrow
