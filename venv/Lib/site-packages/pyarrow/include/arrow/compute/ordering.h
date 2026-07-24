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

#include <string>
#include <vector>

#include "arrow/type.h"
#include "arrow/util/compare.h"
#include "arrow/util/visibility.h"

namespace arrow {
namespace compute {

enum class SortOrder {
  /// Arrange values in increasing order
  Ascending,
  /// Arrange values in decreasing order
  Descending,
};

enum class NullPlacement {
  /// Place nulls and NaNs before any non-null values.
  /// NaNs will come after nulls.
  AtStart,
  /// Place nulls and NaNs after any non-null values.
  /// NaNs will come before nulls.
  AtEnd,
};

/// \brief One sort key for PartitionNthIndices (TODO) and SortIndices
class ARROW_EXPORT SortKey : public util::EqualityComparable<SortKey> {
 public:
  explicit SortKey(FieldRef target, SortOrder order = SortOrder::Ascending)
      : target(std::move(target)), order(order) {}

  bool Equals(const SortKey& other) const;
  std::string ToString() const;

  /// A FieldRef targeting the sort column.
  FieldRef target;
  /// How to order by this sort key.
  SortOrder order;
};

class ARROW_EXPORT Ordering : public util::EqualityComparable<Ordering> {
 public:
  Ordering(std::vector<SortKey> sort_keys,
           NullPlacement null_placement = NullPlacement::AtStart)
      : sort_keys_(std::move(sort_keys)), null_placement_(null_placement) {}
  /// true if data ordered by other is also ordered by this
  ///
  /// For example, if data is ordered by [a, b, c] then it is also ordered
  /// by [a, b] but not by [b, c] or [a, b, c, d].
  ///
  /// [a, b].IsSuborderOf([a, b, c]) - true
  /// [a, b, c].IsSuborderOf([a, b, c]) - true
  /// [b, c].IsSuborderOf([a, b, c]) - false
  /// [a, b, c, d].IsSuborderOf([a, b, c]) - false
  ///
  /// The implicit ordering is not a suborder of any other ordering and
  /// no other ordering is a suborder of it.  The implicit ordering is not a
  /// suborder of itself.
  ///
  /// The unordered ordering is a suborder of all other orderings but no
  /// other ordering is a suborder of it.  The unordered ordering is a suborder
  /// of itself.
  ///
  /// The unordered ordering is a suborder of the implicit ordering.
  bool IsSuborderOf(const Ordering& other) const;

  bool Equals(const Ordering& other) const;
  std::string ToString() const;

  bool is_implicit() const { return is_implicit_; }
  bool is_unordered() const { return !is_implicit_ && sort_keys_.empty(); }

  const std::vector<SortKey>& sort_keys() const { return sort_keys_; }
  NullPlacement null_placement() const { return null_placement_; }

  static const Ordering& Implicit() {
    static const Ordering kImplicit(true);
    return kImplicit;
  }

  static const Ordering& Unordered() {
    static const Ordering kUnordered(false);
    // It is also possible to get an unordered ordering by passing in an empty vector
    // using the normal constructor.  This is ok and useful when ordering comes from user
    // input.
    return kUnordered;
  }

 private:
  explicit Ordering(bool is_implicit)
      : null_placement_(NullPlacement::AtStart), is_implicit_(is_implicit) {}
  /// Column key(s) to order by and how to order by these sort keys.
  std::vector<SortKey> sort_keys_;
  /// Whether nulls and NaNs are placed at the start or at the end
  NullPlacement null_placement_;
  bool is_implicit_ = false;
};

}  // namespace compute
}  // namespace arrow
