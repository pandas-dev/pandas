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
#include <vector>

#include "arrow/compute/kernel.h"
#include "arrow/datum.h"
#include "arrow/result.h"
#include "arrow/util/visibility.h"

namespace arrow {
namespace compute {

/// \brief A segment
/// A segment group is a chunk of continuous rows that have the same segment key. (For
/// example, in ordered time series processing, segment key can be "date", and a segment
/// group can be all the rows that belong to the same date.) A segment group can span
/// across multiple exec batches. A segment is a chunk of continuous rows that has the
/// same segment key within a given batch. When a segment group span cross batches, it
/// will have multiple segments. A segment never spans cross batches. The segment data
/// structure only makes sense when used along with a exec batch.
struct ARROW_EXPORT Segment {
  /// \brief the offset into the batch where the segment starts
  int64_t offset;
  /// \brief the length of the segment
  int64_t length;
  /// \brief whether the segment may be extended by a next one
  bool is_open;
  /// \brief whether the segment extends a preceeding one
  bool extends;
};

inline bool operator==(const Segment& segment1, const Segment& segment2) {
  return segment1.offset == segment2.offset && segment1.length == segment2.length &&
         segment1.is_open == segment2.is_open && segment1.extends == segment2.extends;
}
inline bool operator!=(const Segment& segment1, const Segment& segment2) {
  return !(segment1 == segment2);
}

/// \brief a helper class to divide a batch into segments of equal values
///
/// For example, given a batch with two columns specifed as segment keys:
///
/// A A [other columns]...
/// A A ...
/// A B ...
/// A B ...
/// A A ...
///
/// Then the batch could be divided into 3 segments.  The first would be rows 0 & 1,
/// the second would be rows 2 & 3, and the third would be row 4.
///
/// Further, a segmenter keeps track of the last value seen.  This allows it to calculate
/// segments which span batches.  In our above example the last batch we emit would set
/// the "open" flag, which indicates whether the segment may extend into the next batch.
///
/// If the next call to the segmenter starts with `A A` then that segment would set the
/// "extends" flag, which indicates whether the segment continues the last open batch.
class ARROW_EXPORT RowSegmenter {
 public:
  virtual ~RowSegmenter() = default;

  /// \brief Construct a Segmenter which segments on the specified key types
  ///
  /// \param[in] key_types the specified key types
  /// \param[in] nullable_keys whether values of the specified keys may be null
  /// \param[in] ctx the execution context to use
  static Result<std::unique_ptr<RowSegmenter>> Make(
      const std::vector<TypeHolder>& key_types, bool nullable_keys, ExecContext* ctx);

  /// \brief Return the key types of this segmenter
  virtual const std::vector<TypeHolder>& key_types() const = 0;

  /// \brief Reset this segmenter
  ///
  /// A segmenter normally extends (see `Segment`) a segment from one batch to the next.
  /// If segment-extension is undesirable, for example when each batch is processed
  /// independently, then `Reset` should be invoked before processing the next batch.
  virtual Status Reset() = 0;

  /// \brief Get all segments for the given batch
  virtual Result<std::vector<Segment>> GetSegments(const ExecSpan& batch) = 0;
};

/// Consumes batches of keys and yields batches of the group ids.
class ARROW_EXPORT Grouper {
 public:
  virtual ~Grouper() = default;

  /// Construct a Grouper which receives the specified key types
  static Result<std::unique_ptr<Grouper>> Make(const std::vector<TypeHolder>& key_types,
                                               ExecContext* ctx = default_exec_context());

  /// Reset all intermediate state, make the grouper logically as just `Make`ed.
  /// The underlying buffers, if any, may or may not be released though.
  virtual Status Reset() = 0;

  /// Consume a batch of keys, producing the corresponding group ids as an integer array,
  /// over a slice defined by an offset and length, which defaults to the batch length.
  /// Currently only uint32 indices will be produced, eventually the bit width will only
  /// be as wide as necessary.
  virtual Result<Datum> Consume(const ExecSpan& batch, int64_t offset = 0,
                                int64_t length = -1) = 0;

  /// Like Consume, but groups not already encountered emit null instead of
  /// generating a new group id.
  virtual Result<Datum> Lookup(const ExecSpan& batch, int64_t offset = 0,
                               int64_t length = -1) = 0;

  /// Like Consume, but only populates the Grouper without returning the group ids.
  virtual Status Populate(const ExecSpan& batch, int64_t offset = 0,
                          int64_t length = -1) = 0;

  /// Get current unique keys. May be called multiple times.
  virtual Result<ExecBatch> GetUniques() = 0;

  /// Get the current number of groups.
  virtual uint32_t num_groups() const = 0;

  /// \brief Assemble lists of indices of identical elements.
  ///
  /// \param[in] ids An unsigned, all-valid integral array which will be
  ///                used as grouping criteria.
  /// \param[in] num_groups An upper bound for the elements of ids
  /// \param[in] ctx Execution context to use during the operation
  /// \return A num_groups-long ListArray where the slot at i contains a
  ///         list of indices where i appears in ids.
  ///
  ///   MakeGroupings([
  ///       2,
  ///       2,
  ///       5,
  ///       5,
  ///       2,
  ///       3
  ///   ], 8) == [
  ///       [],
  ///       [],
  ///       [0, 1, 4],
  ///       [5],
  ///       [],
  ///       [2, 3],
  ///       [],
  ///       []
  ///   ]
  static Result<std::shared_ptr<ListArray>> MakeGroupings(
      const UInt32Array& ids, uint32_t num_groups,
      ExecContext* ctx = default_exec_context());

  /// \brief Produce a ListArray whose slots are selections of `array` which correspond to
  /// the provided groupings.
  ///
  /// For example,
  ///   ApplyGroupings([
  ///       [],
  ///       [],
  ///       [0, 1, 4],
  ///       [5],
  ///       [],
  ///       [2, 3],
  ///       [],
  ///       []
  ///   ], [2, 2, 5, 5, 2, 3]) == [
  ///       [],
  ///       [],
  ///       [2, 2, 2],
  ///       [3],
  ///       [],
  ///       [5, 5],
  ///       [],
  ///       []
  ///   ]
  static Result<std::shared_ptr<ListArray>> ApplyGroupings(
      const ListArray& groupings, const Array& array,
      ExecContext* ctx = default_exec_context());
};

}  // namespace compute
}  // namespace arrow
