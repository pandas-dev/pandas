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
#include <utility>

#include "arrow/compute/function_options.h"
#include "arrow/compute/ordering.h"
#include "arrow/result.h"
#include "arrow/type_fwd.h"

namespace arrow {
namespace compute {

class ExecContext;

/// \addtogroup compute-concrete-options
/// @{

class ARROW_EXPORT FilterOptions : public FunctionOptions {
 public:
  /// Configure the action taken when a slot of the selection mask is null
  enum NullSelectionBehavior {
    /// The corresponding filtered value will be removed in the output.
    DROP,
    /// The corresponding filtered value will be null in the output.
    EMIT_NULL,
  };

  explicit FilterOptions(NullSelectionBehavior null_selection = DROP);
  static constexpr char const kTypeName[] = "FilterOptions";
  static FilterOptions Defaults() { return FilterOptions(); }

  NullSelectionBehavior null_selection_behavior = DROP;
};

class ARROW_EXPORT TakeOptions : public FunctionOptions {
 public:
  explicit TakeOptions(bool boundscheck = true);
  static constexpr char const kTypeName[] = "TakeOptions";
  static TakeOptions BoundsCheck() { return TakeOptions(true); }
  static TakeOptions NoBoundsCheck() { return TakeOptions(false); }
  static TakeOptions Defaults() { return BoundsCheck(); }

  bool boundscheck = true;
};

/// \brief Options for the dictionary encode function
class ARROW_EXPORT DictionaryEncodeOptions : public FunctionOptions {
 public:
  /// Configure how null values will be encoded
  enum NullEncodingBehavior {
    /// The null value will be added to the dictionary with a proper index.
    ENCODE,
    /// The null value will be masked in the indices array.
    MASK
  };

  explicit DictionaryEncodeOptions(NullEncodingBehavior null_encoding = MASK);
  static constexpr char const kTypeName[] = "DictionaryEncodeOptions";
  static DictionaryEncodeOptions Defaults() { return DictionaryEncodeOptions(); }

  NullEncodingBehavior null_encoding_behavior = MASK;
};

/// \brief Options for the run-end encode function
class ARROW_EXPORT RunEndEncodeOptions : public FunctionOptions {
 public:
  explicit RunEndEncodeOptions(std::shared_ptr<DataType> run_end_type = int32());
  static constexpr char const kTypeName[] = "RunEndEncodeOptions";
  static RunEndEncodeOptions Defaults() { return RunEndEncodeOptions(); }

  std::shared_ptr<DataType> run_end_type;
};

class ARROW_EXPORT ArraySortOptions : public FunctionOptions {
 public:
  explicit ArraySortOptions(SortOrder order = SortOrder::Ascending,
                            NullPlacement null_placement = NullPlacement::AtEnd);
  static constexpr char const kTypeName[] = "ArraySortOptions";
  static ArraySortOptions Defaults() { return ArraySortOptions(); }

  /// Sorting order
  SortOrder order;
  /// Whether nulls and NaNs are placed at the start or at the end
  NullPlacement null_placement;
};

class ARROW_EXPORT SortOptions : public FunctionOptions {
 public:
  explicit SortOptions(std::vector<SortKey> sort_keys = {},
                       NullPlacement null_placement = NullPlacement::AtEnd);
  explicit SortOptions(const Ordering& ordering);
  static constexpr char const kTypeName[] = "SortOptions";
  static SortOptions Defaults() { return SortOptions(); }
  /// Convenience constructor to create an ordering from SortOptions
  ///
  /// Note: Both classes contain the exact same information.  However,
  /// sort_options should only be used in a "function options" context while Ordering
  /// is used more generally.
  Ordering AsOrdering() && { return Ordering(std::move(sort_keys), null_placement); }
  Ordering AsOrdering() const& { return Ordering(sort_keys, null_placement); }

  /// Column key(s) to order by and how to order by these sort keys.
  std::vector<SortKey> sort_keys;
  /// Whether nulls and NaNs are placed at the start or at the end
  NullPlacement null_placement;
};

/// \brief SelectK options
class ARROW_EXPORT SelectKOptions : public FunctionOptions {
 public:
  explicit SelectKOptions(int64_t k = -1, std::vector<SortKey> sort_keys = {});
  static constexpr char const kTypeName[] = "SelectKOptions";
  static SelectKOptions Defaults() { return SelectKOptions(); }

  static SelectKOptions TopKDefault(int64_t k, std::vector<std::string> key_names = {}) {
    std::vector<SortKey> keys;
    for (const auto& name : key_names) {
      keys.emplace_back(SortKey(name, SortOrder::Descending));
    }
    if (key_names.empty()) {
      keys.emplace_back(SortKey("not-used", SortOrder::Descending));
    }
    return SelectKOptions{k, keys};
  }
  static SelectKOptions BottomKDefault(int64_t k,
                                       std::vector<std::string> key_names = {}) {
    std::vector<SortKey> keys;
    for (const auto& name : key_names) {
      keys.emplace_back(SortKey(name, SortOrder::Ascending));
    }
    if (key_names.empty()) {
      keys.emplace_back(SortKey("not-used", SortOrder::Ascending));
    }
    return SelectKOptions{k, keys};
  }

  /// The number of `k` elements to keep.
  int64_t k;
  /// Column key(s) to order by and how to order by these sort keys.
  std::vector<SortKey> sort_keys;
};

/// \brief Rank options
class ARROW_EXPORT RankOptions : public FunctionOptions {
 public:
  /// Configure how ties between equal values are handled
  enum Tiebreaker {
    /// Ties get the smallest possible rank in sorted order.
    Min,
    /// Ties get the largest possible rank in sorted order.
    Max,
    /// Ranks are assigned in order of when ties appear in the input.
    /// This ensures the ranks are a stable permutation of the input.
    First,
    /// The ranks span a dense [1, M] interval where M is the number
    /// of distinct values in the input.
    Dense
  };

  explicit RankOptions(std::vector<SortKey> sort_keys = {},
                       NullPlacement null_placement = NullPlacement::AtEnd,
                       Tiebreaker tiebreaker = RankOptions::First);
  /// Convenience constructor for array inputs
  explicit RankOptions(SortOrder order,
                       NullPlacement null_placement = NullPlacement::AtEnd,
                       Tiebreaker tiebreaker = RankOptions::First)
      : RankOptions({SortKey("", order)}, null_placement, tiebreaker) {}

  static constexpr char const kTypeName[] = "RankOptions";
  static RankOptions Defaults() { return RankOptions(); }

  /// Column key(s) to order by and how to order by these sort keys.
  std::vector<SortKey> sort_keys;
  /// Whether nulls and NaNs are placed at the start or at the end
  NullPlacement null_placement;
  /// Tiebreaker for dealing with equal values in ranks
  Tiebreaker tiebreaker;
};

/// \brief Quantile rank options
class ARROW_EXPORT RankQuantileOptions : public FunctionOptions {
 public:
  explicit RankQuantileOptions(std::vector<SortKey> sort_keys = {},
                               NullPlacement null_placement = NullPlacement::AtEnd);
  /// Convenience constructor for array inputs
  explicit RankQuantileOptions(SortOrder order,
                               NullPlacement null_placement = NullPlacement::AtEnd)
      : RankQuantileOptions({SortKey("", order)}, null_placement) {}

  static constexpr char const kTypeName[] = "RankQuantileOptions";
  static RankQuantileOptions Defaults() { return RankQuantileOptions(); }

  /// Column key(s) to order by and how to order by these sort keys.
  std::vector<SortKey> sort_keys;
  /// Whether nulls and NaNs are placed at the start or at the end
  NullPlacement null_placement;
};

/// \brief Partitioning options for NthToIndices
class ARROW_EXPORT PartitionNthOptions : public FunctionOptions {
 public:
  explicit PartitionNthOptions(int64_t pivot,
                               NullPlacement null_placement = NullPlacement::AtEnd);
  PartitionNthOptions() : PartitionNthOptions(0) {}
  static constexpr char const kTypeName[] = "PartitionNthOptions";

  /// The index into the equivalent sorted array of the partition pivot element.
  int64_t pivot;
  /// Whether nulls and NaNs are partitioned at the start or at the end
  NullPlacement null_placement;
};

class ARROW_EXPORT WinsorizeOptions : public FunctionOptions {
 public:
  WinsorizeOptions(double lower_limit, double upper_limit);
  WinsorizeOptions() : WinsorizeOptions(0, 1) {}
  static constexpr char const kTypeName[] = "WinsorizeOptions";

  /// The quantile below which all values are replaced with the quantile's value.
  ///
  /// For example, if lower_limit = 0.05, then all values in the lower 5% percentile
  /// will be replaced with the 5% percentile value.
  double lower_limit;

  /// The quantile above which all values are replaced with the quantile's value.
  ///
  /// For example, if upper_limit = 0.95, then all values in the upper 95% percentile
  /// will be replaced with the 95% percentile value.
  double upper_limit;
};

/// \brief Options for cumulative functions
/// \note Also aliased as CumulativeSumOptions for backward compatibility
class ARROW_EXPORT CumulativeOptions : public FunctionOptions {
 public:
  explicit CumulativeOptions(bool skip_nulls = false);
  explicit CumulativeOptions(double start, bool skip_nulls = false);
  explicit CumulativeOptions(std::shared_ptr<Scalar> start, bool skip_nulls = false);
  static constexpr char const kTypeName[] = "CumulativeOptions";
  static CumulativeOptions Defaults() { return CumulativeOptions(); }

  /// Optional starting value for cumulative operation computation, default depends on the
  /// operation and input type.
  /// - sum: 0
  /// - prod: 1
  /// - min: maximum of the input type
  /// - max: minimum of the input type
  /// - mean: start is ignored because it has no meaning for mean
  std::optional<std::shared_ptr<Scalar>> start;

  /// If true, nulls in the input are ignored and produce a corresponding null output.
  /// When false, the first null encountered is propagated through the remaining output.
  bool skip_nulls = false;
};
using CumulativeSumOptions = CumulativeOptions;  // For backward compatibility

/// \brief Options for pairwise functions
class ARROW_EXPORT PairwiseOptions : public FunctionOptions {
 public:
  explicit PairwiseOptions(int64_t periods = 1);
  static constexpr char const kTypeName[] = "PairwiseOptions";
  static PairwiseOptions Defaults() { return PairwiseOptions(); }

  /// Periods to shift for applying the binary operation, accepts negative values.
  int64_t periods = 1;
};

/// \brief Options for list_flatten function
class ARROW_EXPORT ListFlattenOptions : public FunctionOptions {
 public:
  explicit ListFlattenOptions(bool recursive = false);
  static constexpr char const kTypeName[] = "ListFlattenOptions";
  static ListFlattenOptions Defaults() { return ListFlattenOptions(); }

  /// \brief If true, the list is flattened recursively until a non-list
  /// array is formed.
  bool recursive = false;
};

/// \brief Options for inverse_permutation function
class ARROW_EXPORT InversePermutationOptions : public FunctionOptions {
 public:
  explicit InversePermutationOptions(int64_t max_index = -1,
                                     std::shared_ptr<DataType> output_type = NULLPTR);
  static constexpr char const kTypeName[] = "InversePermutationOptions";
  static InversePermutationOptions Defaults() { return InversePermutationOptions(); }

  /// \brief The max value in the input indices to allow. The length of the function's
  /// output will be this value plus 1. If negative, this value will be set to the length
  /// of the input indices minus 1 and the length of the function's output will be the
  /// length of the input indices.
  int64_t max_index = -1;
  /// \brief The type of the output inverse permutation. If null, the output will be of
  /// the same type as the input indices, otherwise must be signed integer type. An
  /// invalid error will be reported if this type is not able to store the length of the
  /// input indices.
  std::shared_ptr<DataType> output_type = NULLPTR;
};

/// \brief Options for scatter function
class ARROW_EXPORT ScatterOptions : public FunctionOptions {
 public:
  explicit ScatterOptions(int64_t max_index = -1);
  static constexpr char const kTypeName[] = "ScatterOptions";
  static ScatterOptions Defaults() { return ScatterOptions(); }

  /// \brief The max value in the input indices to allow. The length of the function's
  /// output will be this value plus 1. If negative, this value will be set to the length
  /// of the input indices minus 1 and the length of the function's output will be the
  /// length of the input indices.
  int64_t max_index = -1;
};

/// @}

/// \brief Filter with a boolean selection filter
///
/// The output will be populated with values from the input at positions
/// where the selection filter is not 0. Nulls in the filter will be handled
/// based on options.null_selection_behavior.
///
/// For example given values = ["a", "b", "c", null, "e", "f"] and
/// filter = [0, 1, 1, 0, null, 1], the output will be
/// (null_selection_behavior == DROP)      = ["b", "c", "f"]
/// (null_selection_behavior == EMIT_NULL) = ["b", "c", null, "f"]
///
/// \param[in] values array to filter
/// \param[in] filter indicates which values should be filtered out
/// \param[in] options configures null_selection_behavior
/// \param[in] ctx the function execution context, optional
/// \return the resulting datum
ARROW_EXPORT
Result<Datum> Filter(const Datum& values, const Datum& filter,
                     const FilterOptions& options = FilterOptions::Defaults(),
                     ExecContext* ctx = NULLPTR);

namespace internal {

// These internal functions are implemented in kernels/vector_selection.cc

/// \brief Return the number of selected indices in the boolean filter
///
/// \param filter a plain or run-end encoded boolean array with or without nulls
/// \param null_selection how to handle nulls in the filter
ARROW_EXPORT
int64_t GetFilterOutputSize(const ArraySpan& filter,
                            FilterOptions::NullSelectionBehavior null_selection);

/// \brief Compute uint64 selection indices for use with Take given a boolean
/// filter
///
/// \param filter a plain or run-end encoded boolean array with or without nulls
/// \param null_selection how to handle nulls in the filter
ARROW_EXPORT
Result<std::shared_ptr<ArrayData>> GetTakeIndices(
    const ArraySpan& filter, FilterOptions::NullSelectionBehavior null_selection,
    MemoryPool* memory_pool = default_memory_pool());

}  // namespace internal

/// \brief ReplaceWithMask replaces each value in the array corresponding
/// to a true value in the mask with the next element from `replacements`.
///
/// \param[in] values Array input to replace
/// \param[in] mask Array or Scalar of Boolean mask values
/// \param[in] replacements The replacement values to draw from. There must
/// be as many replacement values as true values in the mask.
/// \param[in] ctx the function execution context, optional
///
/// \return the resulting datum
///
/// \since 5.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> ReplaceWithMask(const Datum& values, const Datum& mask,
                              const Datum& replacements, ExecContext* ctx = NULLPTR);

/// \brief FillNullForward fill null values in forward direction
///
/// The output array will be of the same type as the input values
/// array, with replaced null values in forward direction.
///
/// For example given values = ["a", "b", "c", null, null, "f"],
/// the output will be = ["a", "b", "c", "c", "c", "f"]
///
/// \param[in] values datum from which to take
/// \param[in] ctx the function execution context, optional
/// \return the resulting datum
ARROW_EXPORT
Result<Datum> FillNullForward(const Datum& values, ExecContext* ctx = NULLPTR);

/// \brief FillNullBackward fill null values in backward direction
///
/// The output array will be of the same type as the input values
/// array, with replaced null values in backward direction.
///
/// For example given values = ["a", "b", "c", null, null, "f"],
/// the output will be = ["a", "b", "c", "f", "f", "f"]
///
/// \param[in] values datum from which to take
/// \param[in] ctx the function execution context, optional
/// \return the resulting datum
ARROW_EXPORT
Result<Datum> FillNullBackward(const Datum& values, ExecContext* ctx = NULLPTR);

/// \brief Take from an array of values at indices in another array
///
/// The output array will be of the same type as the input values
/// array, with elements taken from the values array at the given
/// indices. If an index is null then the taken element will be null.
///
/// For example given values = ["a", "b", "c", null, "e", "f"] and
/// indices = [2, 1, null, 3], the output will be
/// = [values[2], values[1], null, values[3]]
/// = ["c", "b", null, null]
///
/// \param[in] values datum from which to take
/// \param[in] indices which values to take
/// \param[in] options options
/// \param[in] ctx the function execution context, optional
/// \return the resulting datum
ARROW_EXPORT
Result<Datum> Take(const Datum& values, const Datum& indices,
                   const TakeOptions& options = TakeOptions::Defaults(),
                   ExecContext* ctx = NULLPTR);

/// \brief Take with Array inputs and output
ARROW_EXPORT
Result<std::shared_ptr<Array>> Take(const Array& values, const Array& indices,
                                    const TakeOptions& options = TakeOptions::Defaults(),
                                    ExecContext* ctx = NULLPTR);

/// \brief Drop Null from an array of values
///
/// The output array will be of the same type as the input values
/// array, with elements taken from the values array without nulls.
///
/// For example given values = ["a", "b", "c", null, "e", "f"],
/// the output will be = ["a", "b", "c", "e", "f"]
///
/// \param[in] values datum from which to take
/// \param[in] ctx the function execution context, optional
/// \return the resulting datum
ARROW_EXPORT
Result<Datum> DropNull(const Datum& values, ExecContext* ctx = NULLPTR);

/// \brief DropNull with Array inputs and output
ARROW_EXPORT
Result<std::shared_ptr<Array>> DropNull(const Array& values, ExecContext* ctx = NULLPTR);

/// \brief Return indices that partition an array around n-th sorted element.
///
/// Find index of n-th(0 based) smallest value and perform indirect
/// partition of an array around that element. Output indices[0 ~ n-1]
/// holds values no greater than n-th element, and indices[n+1 ~ end]
/// holds values no less than n-th element. Elements in each partition
/// is not sorted. Nulls will be partitioned to the end of the output.
/// Output is not guaranteed to be stable.
///
/// \param[in] values array to be partitioned
/// \param[in] n pivot array around sorted n-th element
/// \param[in] ctx the function execution context, optional
/// \return offsets indices that would partition an array
ARROW_EXPORT
Result<std::shared_ptr<Array>> NthToIndices(const Array& values, int64_t n,
                                            ExecContext* ctx = NULLPTR);

/// \brief Return indices that partition an array around n-th sorted element.
///
/// This overload takes a PartitionNthOptions specifying the pivot index
/// and the null handling.
///
/// \param[in] values array to be partitioned
/// \param[in] options options including pivot index and null handling
/// \param[in] ctx the function execution context, optional
/// \return offsets indices that would partition an array
ARROW_EXPORT
Result<std::shared_ptr<Array>> NthToIndices(const Array& values,
                                            const PartitionNthOptions& options,
                                            ExecContext* ctx = NULLPTR);

/// \brief Return indices that would select the first `k` elements.
///
/// Perform an indirect sort of the datum, keeping only the first `k` elements. The output
/// array will contain indices such that the item indicated by the k-th index will be in
/// the position it would be if the datum were sorted by `options.sort_keys`. However,
/// indices of null values will not be part of the output. The sort is not guaranteed to
/// be stable.
///
/// \param[in] datum datum to be partitioned
/// \param[in] options options
/// \param[in] ctx the function execution context, optional
/// \return a datum with the same schema as the input
ARROW_EXPORT
Result<std::shared_ptr<Array>> SelectKUnstable(const Datum& datum,
                                               const SelectKOptions& options,
                                               ExecContext* ctx = NULLPTR);

/// \brief Return the indices that would sort an array.
///
/// Perform an indirect sort of array. The output array will contain
/// indices that would sort an array, which would be the same length
/// as input. Nulls will be stably partitioned to the end of the output
/// regardless of order.
///
/// For example given array = [null, 1, 3.3, null, 2, 5.3] and order
/// = SortOrder::DESCENDING, the output will be [5, 2, 4, 1, 0,
/// 3].
///
/// \param[in] array array to sort
/// \param[in] order ascending or descending
/// \param[in] ctx the function execution context, optional
/// \return offsets indices that would sort an array
ARROW_EXPORT
Result<std::shared_ptr<Array>> SortIndices(const Array& array,
                                           SortOrder order = SortOrder::Ascending,
                                           ExecContext* ctx = NULLPTR);

/// \brief Return the indices that would sort an array.
///
/// This overload takes a ArraySortOptions specifying the sort order
/// and the null handling.
///
/// \param[in] array array to sort
/// \param[in] options options including sort order and null handling
/// \param[in] ctx the function execution context, optional
/// \return offsets indices that would sort an array
ARROW_EXPORT
Result<std::shared_ptr<Array>> SortIndices(const Array& array,
                                           const ArraySortOptions& options,
                                           ExecContext* ctx = NULLPTR);

/// \brief Return the indices that would sort a chunked array.
///
/// Perform an indirect sort of chunked array. The output array will
/// contain indices that would sort a chunked array, which would be
/// the same length as input. Nulls will be stably partitioned to the
/// end of the output regardless of order.
///
/// For example given chunked_array = [[null, 1], [3.3], [null, 2,
/// 5.3]] and order = SortOrder::DESCENDING, the output will be [5, 2,
/// 4, 1, 0, 3].
///
/// \param[in] chunked_array chunked array to sort
/// \param[in] order ascending or descending
/// \param[in] ctx the function execution context, optional
/// \return offsets indices that would sort an array
ARROW_EXPORT
Result<std::shared_ptr<Array>> SortIndices(const ChunkedArray& chunked_array,
                                           SortOrder order = SortOrder::Ascending,
                                           ExecContext* ctx = NULLPTR);

/// \brief Return the indices that would sort a chunked array.
///
/// This overload takes a ArraySortOptions specifying the sort order
/// and the null handling.
///
/// \param[in] chunked_array chunked array to sort
/// \param[in] options options including sort order and null handling
/// \param[in] ctx the function execution context, optional
/// \return offsets indices that would sort an array
ARROW_EXPORT
Result<std::shared_ptr<Array>> SortIndices(const ChunkedArray& chunked_array,
                                           const ArraySortOptions& options,
                                           ExecContext* ctx = NULLPTR);

/// \brief Return the indices that would sort an input in the
/// specified order. Input is one of array, chunked array record batch
/// or table.
///
/// Perform an indirect sort of input. The output array will contain
/// indices that would sort an input, which would be the same length
/// as input. Nulls will be stably partitioned to the start or to the end
/// of the output depending on SortOrder::null_placement.
///
/// For example given input (table) = {
/// "column1": [[null,   1], [   3, null, 2, 1]],
/// "column2": [[   5], [3,   null, null, 5, 5]],
/// } and options = {
/// {"column1", SortOrder::Ascending},
/// {"column2", SortOrder::Descending},
/// }, the output will be [5, 1, 4, 2, 0, 3].
///
/// \param[in] datum array, chunked array, record batch or table to sort
/// \param[in] options options
/// \param[in] ctx the function execution context, optional
/// \return offsets indices that would sort a table
ARROW_EXPORT
Result<std::shared_ptr<Array>> SortIndices(const Datum& datum, const SortOptions& options,
                                           ExecContext* ctx = NULLPTR);

/// \brief Compute unique elements from an array-like object
///
/// Note if a null occurs in the input it will NOT be included in the output.
///
/// \param[in] datum array-like input
/// \param[in] ctx the function execution context, optional
/// \return result as Array
///
/// \since 1.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<std::shared_ptr<Array>> Unique(const Datum& datum, ExecContext* ctx = NULLPTR);

// Constants for accessing the output of ValueCounts
ARROW_EXPORT extern const char kValuesFieldName[];
ARROW_EXPORT extern const char kCountsFieldName[];
ARROW_EXPORT extern const int32_t kValuesFieldIndex;
ARROW_EXPORT extern const int32_t kCountsFieldIndex;

/// \brief Return counts of unique elements from an array-like object.
///
/// Note that the counts do not include counts for nulls in the array.  These can be
/// obtained separately from metadata.
///
/// For floating point arrays there is no attempt to normalize -0.0, 0.0 and NaN values
/// which can lead to unexpected results if the input Array has these values.
///
/// \param[in] value array-like input
/// \param[in] ctx the function execution context, optional
/// \return counts An array of  <input type "Values", int64_t "Counts"> structs.
///
/// \since 1.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<std::shared_ptr<StructArray>> ValueCounts(const Datum& value,
                                                 ExecContext* ctx = NULLPTR);

/// \brief Dictionary-encode values in an array-like object
///
/// Any nulls encountered in the dictionary will be handled according to the
/// specified null encoding behavior.
///
/// For example, given values ["a", "b", null, "a", null] the output will be
/// (null_encoding == ENCODE) Indices: [0, 1, 2, 0, 2] / Dict: ["a", "b", null]
/// (null_encoding == MASK)   Indices: [0, 1, null, 0, null] / Dict: ["a", "b"]
///
/// If the input is already dictionary encoded this function is a no-op unless
/// it needs to modify the null_encoding (TODO)
///
/// \param[in] data array-like input
/// \param[in] ctx the function execution context, optional
/// \param[in] options configures null encoding behavior
/// \return result with same shape and type as input
///
/// \since 1.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> DictionaryEncode(
    const Datum& data,
    const DictionaryEncodeOptions& options = DictionaryEncodeOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Run-end-encode values in an array-like object
///
/// The returned run-end encoded type uses the same value type of the input and
/// run-end type defined in the options.
///
/// \param[in] value array-like input
/// \param[in] options configures encoding behavior
/// \param[in] ctx the function execution context, optional
/// \return result with same shape but run-end encoded
///
/// \since 12.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> RunEndEncode(
    const Datum& value,
    const RunEndEncodeOptions& options = RunEndEncodeOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Decode a Run-End Encoded array to a plain array
///
/// The output data type is the same as the values array type of run-end encoded
/// input.
///
/// \param[in] value run-end-encoded input
/// \param[in] ctx the function execution context, optional
/// \return plain array resulting from decoding the run-end encoded input
///
/// \since 12.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> RunEndDecode(const Datum& value, ExecContext* ctx = NULLPTR);

/// \brief Compute the cumulative sum of an array-like object
///
/// \param[in] values array-like input
/// \param[in] options configures cumulative sum behavior
/// \param[in] check_overflow whether to check for overflow, if true, return Invalid
/// status on overflow, otherwise wrap around on overflow
/// \param[in] ctx the function execution context, optional
ARROW_EXPORT
Result<Datum> CumulativeSum(
    const Datum& values, const CumulativeOptions& options = CumulativeOptions::Defaults(),
    bool check_overflow = false, ExecContext* ctx = NULLPTR);

/// \brief Compute the cumulative product of an array-like object
///
/// \param[in] values array-like input
/// \param[in] options configures cumulative prod behavior
/// \param[in] check_overflow whether to check for overflow, if true, return Invalid
/// status on overflow, otherwise wrap around on overflow
/// \param[in] ctx the function execution context, optional
ARROW_EXPORT
Result<Datum> CumulativeProd(
    const Datum& values, const CumulativeOptions& options = CumulativeOptions::Defaults(),
    bool check_overflow = false, ExecContext* ctx = NULLPTR);

/// \brief Compute the cumulative max of an array-like object
///
/// \param[in] values array-like input
/// \param[in] options configures cumulative max behavior
/// \param[in] ctx the function execution context, optional
ARROW_EXPORT
Result<Datum> CumulativeMax(
    const Datum& values, const CumulativeOptions& options = CumulativeOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Compute the cumulative min of an array-like object
///
/// \param[in] values array-like input
/// \param[in] options configures cumulative min behavior
/// \param[in] ctx the function execution context, optional
ARROW_EXPORT
Result<Datum> CumulativeMin(
    const Datum& values, const CumulativeOptions& options = CumulativeOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Compute the cumulative mean of an array-like object
///
/// \param[in] values array-like input
/// \param[in] options configures cumulative mean behavior, `start` is ignored
/// \param[in] ctx the function execution context, optional
ARROW_EXPORT
Result<Datum> CumulativeMean(
    const Datum& values, const CumulativeOptions& options = CumulativeOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Return the first order difference of an array.
///
/// Computes the first order difference of an array, i.e.
///   output[i] = input[i] - input[i - p]  if i >= p
///   output[i] = null                     otherwise
/// where p is the period. For example, with p = 1,
///   Diff([1, 4, 9, 10, 15]) = [null, 3, 5, 1, 5].
/// With p = 2,
///   Diff([1, 4, 9, 10, 15]) = [null, null, 8, 6, 6]
/// p can also be negative, in which case the diff is computed in
/// the opposite direction.
/// \param[in] array array input
/// \param[in] options options, specifying overflow behavior and period
/// \param[in] check_overflow whether to return error on overflow
/// \param[in] ctx the function execution context, optional
/// \return result as array
ARROW_EXPORT
Result<std::shared_ptr<Array>> PairwiseDiff(const Array& array,
                                            const PairwiseOptions& options,
                                            bool check_overflow = false,
                                            ExecContext* ctx = NULLPTR);

/// \brief Return the inverse permutation of the given indices.
///
/// For indices[i] = x, inverse_permutation[x] = i. And inverse_permutation[x] = null if x
/// does not appear in the input indices. Indices must be in the range of [0, max_index],
/// or null, which will be ignored. If multiple indices point to the same value, the last
/// one is used.
///
/// For example, with
///   indices = [null, 0, null, 2, 4, 1, 1]
/// the inverse permutation is
///   [1, 6, 3, null, 4, null, null]
/// if max_index = 6.
///
/// \param[in] indices array-like indices
/// \param[in] options configures the max index and the output type
/// \param[in] ctx the function execution context, optional
/// \return the resulting inverse permutation
///
/// \since 20.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> InversePermutation(
    const Datum& indices,
    const InversePermutationOptions& options = InversePermutationOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Scatter the values into specified positions according to the indices.
///
/// For indices[i] = x, output[x] = values[i]. And output[x] = null if x does not appear
/// in the input indices. Indices must be in the range of [0, max_index], or null, in
/// which case the corresponding value will be ignored. If multiple indices point to the
/// same value, the last one is used.
///
/// For example, with
///   values = [a, b, c, d, e, f, g]
///   indices = [null, 0, null, 2, 4, 1, 1]
/// the output is
///   [b, g, d, null, e, null, null]
/// if max_index = 6.
///
/// \param[in] values datum to scatter
/// \param[in] indices array-like indices
/// \param[in] options configures the max index of to scatter
/// \param[in] ctx the function execution context, optional
/// \return the resulting datum
///
/// \since 20.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Scatter(const Datum& values, const Datum& indices,
                      const ScatterOptions& options = ScatterOptions::Defaults(),
                      ExecContext* ctx = NULLPTR);

}  // namespace compute
}  // namespace arrow
