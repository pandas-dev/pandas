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

// Eager evaluation convenience APIs for invoking common functions, including
// necessary memory allocations

#pragma once

#include <vector>

#include "arrow/compute/function_options.h"
#include "arrow/datum.h"
#include "arrow/result.h"
#include "arrow/util/macros.h"
#include "arrow/util/visibility.h"

namespace arrow {

class Array;

namespace compute {

class ExecContext;

// ----------------------------------------------------------------------
// Aggregate functions

/// \addtogroup compute-concrete-options
/// @{

/// \brief Control general scalar aggregate kernel behavior
///
/// By default, null values are ignored (skip_nulls = true).
class ARROW_EXPORT ScalarAggregateOptions : public FunctionOptions {
 public:
  explicit ScalarAggregateOptions(bool skip_nulls = true, uint32_t min_count = 1);
  static constexpr char const kTypeName[] = "ScalarAggregateOptions";
  static ScalarAggregateOptions Defaults() { return ScalarAggregateOptions{}; }

  /// If true (the default), null values are ignored. Otherwise, if any value is null,
  /// emit null.
  bool skip_nulls;
  /// If less than this many non-null values are observed, emit null.
  uint32_t min_count;
};

/// \brief Control count aggregate kernel behavior.
///
/// By default, only non-null values are counted.
class ARROW_EXPORT CountOptions : public FunctionOptions {
 public:
  enum CountMode {
    /// Count only non-null values.
    ONLY_VALID = 0,
    /// Count only null values.
    ONLY_NULL,
    /// Count both non-null and null values.
    ALL,
  };
  explicit CountOptions(CountMode mode = CountMode::ONLY_VALID);
  static constexpr char const kTypeName[] = "CountOptions";
  static CountOptions Defaults() { return CountOptions{}; }

  CountMode mode;
};

/// \brief Control Mode kernel behavior
///
/// Returns top-n common values and counts.
/// By default, returns the most common value and count.
class ARROW_EXPORT ModeOptions : public FunctionOptions {
 public:
  explicit ModeOptions(int64_t n = 1, bool skip_nulls = true, uint32_t min_count = 0);
  static constexpr char const kTypeName[] = "ModeOptions";
  static ModeOptions Defaults() { return ModeOptions{}; }

  int64_t n = 1;
  /// If true (the default), null values are ignored. Otherwise, if any value is null,
  /// emit null.
  bool skip_nulls;
  /// If less than this many non-null values are observed, emit null.
  uint32_t min_count;
};

/// \brief Control Delta Degrees of Freedom (ddof) of Variance and Stddev kernel
///
/// The divisor used in calculations is N - ddof, where N is the number of elements.
/// By default, ddof is zero, and population variance or stddev is returned.
class ARROW_EXPORT VarianceOptions : public FunctionOptions {
 public:
  explicit VarianceOptions(int ddof = 0, bool skip_nulls = true, uint32_t min_count = 0);
  static constexpr char const kTypeName[] = "VarianceOptions";
  static VarianceOptions Defaults() { return VarianceOptions{}; }

  int ddof = 0;
  /// If true (the default), null values are ignored. Otherwise, if any value is null,
  /// emit null.
  bool skip_nulls;
  /// If less than this many non-null values are observed, emit null.
  uint32_t min_count;
};

/// \brief Control Skew and Kurtosis kernel behavior
class ARROW_EXPORT SkewOptions : public FunctionOptions {
 public:
  explicit SkewOptions(bool skip_nulls = true, bool biased = true,
                       uint32_t min_count = 0);
  static constexpr char const kTypeName[] = "SkewOptions";
  static SkewOptions Defaults() { return SkewOptions{}; }

  /// If true (the default), null values are ignored. Otherwise, if any value is null,
  /// emit null.
  bool skip_nulls;
  /// If true (the default), the calculated value is biased. If false, the calculated
  /// value includes a correction factor to reduce bias, making it more accurate for
  /// small sample sizes.
  bool biased;
  /// If less than this many non-null values are observed, emit null.
  uint32_t min_count;
};

/// \brief Control Quantile kernel behavior
///
/// By default, returns the median value.
class ARROW_EXPORT QuantileOptions : public FunctionOptions {
 public:
  /// Interpolation method to use when quantile lies between two data points
  enum Interpolation {
    LINEAR = 0,
    LOWER,
    HIGHER,
    NEAREST,
    MIDPOINT,
  };

  explicit QuantileOptions(double q = 0.5, enum Interpolation interpolation = LINEAR,
                           bool skip_nulls = true, uint32_t min_count = 0);

  explicit QuantileOptions(std::vector<double> q,
                           enum Interpolation interpolation = LINEAR,
                           bool skip_nulls = true, uint32_t min_count = 0);

  static constexpr char const kTypeName[] = "QuantileOptions";
  static QuantileOptions Defaults() { return QuantileOptions{}; }

  /// probability level of quantile must be between 0 and 1 inclusive
  std::vector<double> q;
  enum Interpolation interpolation;
  /// If true (the default), null values are ignored. Otherwise, if any value is null,
  /// emit null.
  bool skip_nulls;
  /// If less than this many non-null values are observed, emit null.
  uint32_t min_count;
};

/// \brief Control TDigest approximate quantile kernel behavior
///
/// By default, returns the median value.
class ARROW_EXPORT TDigestOptions : public FunctionOptions {
 public:
  explicit TDigestOptions(double q = 0.5, uint32_t delta = 100,
                          uint32_t buffer_size = 500, bool skip_nulls = true,
                          uint32_t min_count = 0);
  explicit TDigestOptions(std::vector<double> q, uint32_t delta = 100,
                          uint32_t buffer_size = 500, bool skip_nulls = true,
                          uint32_t min_count = 0);
  static constexpr char const kTypeName[] = "TDigestOptions";
  static TDigestOptions Defaults() { return TDigestOptions{}; }

  /// probability level of quantile must be between 0 and 1 inclusive
  std::vector<double> q;
  /// compression parameter, default 100
  uint32_t delta;
  /// input buffer size, default 500
  uint32_t buffer_size;
  /// If true (the default), null values are ignored. Otherwise, if any value is null,
  /// emit null.
  bool skip_nulls;
  /// If less than this many non-null values are observed, emit null.
  uint32_t min_count;
};

/// \brief Control Pivot kernel behavior
///
/// These options apply to the "pivot_wider" and "hash_pivot_wider" functions.
///
/// Constraints:
/// - The corresponding `Aggregate::target` must have two FieldRef elements;
///   the first one points to the pivot key column, the second points to the
///   pivoted data column.
/// - The pivot key column can be string, binary or integer; its values will be
///   matched against `key_names` in order to dispatch the pivoted data into
///   the output. If the pivot key column is not string-like, the `key_names`
///   will be cast to the pivot key type.
///
/// "pivot_wider" example
/// ---------------------
///
/// Assuming the following two input columns with types utf8 and int16 (respectively):
/// ```
/// width   |  11
/// height  |  13
/// ```
/// and the options `PivotWiderOptions(.key_names = {"height", "width"})`
///
/// then the output will be a scalar with the type
/// `struct{"height": int16, "width": int16}`
/// and the value `{"height": 13, "width": 11}`.
///
/// "hash_pivot_wider" example
/// --------------------------
///
/// Assuming the following input with schema
/// `{"group": int32, "key": utf8, "value": int16}`:
/// ```
///  group |  key     |  value
/// -----------------------------
///   1    |  height  |    11
///   1    |  width   |    12
///   2    |  width   |    13
///   3    |  height  |    14
///   3    |  depth   |    15
/// ```
/// and the following settings:
/// - a hash grouping key "group"
/// - Aggregate(
///     .function = "hash_pivot_wider",
///     .options = PivotWiderOptions(.key_names = {"height", "width"}),
///     .target = {"key", "value"},
///     .name = {"properties"})
///
/// then the output will have the schema
/// `{"group": int32, "properties": struct{"height": int16, "width": int16}}`
/// and the following value:
/// ```
///  group |     properties
///        |  height  |   width
/// -----------------------------
///   1    |   11     |    12
///   2    |   null   |    13
///   3    |   14     |    null
/// ```
class ARROW_EXPORT PivotWiderOptions : public FunctionOptions {
 public:
  /// Configure the behavior of pivot keys not in `key_names`
  enum UnexpectedKeyBehavior {
    /// Unexpected pivot keys are ignored silently
    kIgnore,
    /// Unexpected pivot keys return a KeyError
    kRaise
  };

  explicit PivotWiderOptions(std::vector<std::string> key_names,
                             UnexpectedKeyBehavior unexpected_key_behavior = kIgnore);
  // Default constructor for serialization
  PivotWiderOptions();
  static constexpr char const kTypeName[] = "PivotWiderOptions";
  static PivotWiderOptions Defaults() { return PivotWiderOptions{}; }

  /// The values expected in the pivot key column
  std::vector<std::string> key_names;
  /// The behavior when pivot keys not in `key_names` are encountered
  UnexpectedKeyBehavior unexpected_key_behavior = kIgnore;
};

/// \brief Control Index kernel behavior
class ARROW_EXPORT IndexOptions : public FunctionOptions {
 public:
  explicit IndexOptions(std::shared_ptr<Scalar> value);
  // Default constructor for serialization
  IndexOptions();
  static constexpr char const kTypeName[] = "IndexOptions";

  std::shared_ptr<Scalar> value;
};

/// \brief Configure a grouped aggregation
struct ARROW_EXPORT Aggregate {
  Aggregate() = default;

  Aggregate(std::string function, std::shared_ptr<FunctionOptions> options,
            std::vector<FieldRef> target, std::string name = "")
      : function(std::move(function)),
        options(std::move(options)),
        target(std::move(target)),
        name(std::move(name)) {}

  Aggregate(std::string function, std::shared_ptr<FunctionOptions> options,
            FieldRef target, std::string name = "")
      : Aggregate(std::move(function), std::move(options),
                  std::vector<FieldRef>{std::move(target)}, std::move(name)) {}

  Aggregate(std::string function, FieldRef target, std::string name)
      : Aggregate(std::move(function), /*options=*/NULLPTR,
                  std::vector<FieldRef>{std::move(target)}, std::move(name)) {}

  Aggregate(std::string function, std::string name)
      : Aggregate(std::move(function), /*options=*/NULLPTR,
                  /*target=*/std::vector<FieldRef>{}, std::move(name)) {}

  /// the name of the aggregation function
  std::string function;

  /// options for the aggregation function
  std::shared_ptr<FunctionOptions> options;

  /// zero or more fields to which aggregations will be applied
  std::vector<FieldRef> target;

  /// optional output field name for aggregations
  std::string name;
};

/// @}

/// \brief Count values in an array.
///
/// \param[in] options counting options, see CountOptions for more information
/// \param[in] datum to count
/// \param[in] ctx the function execution context, optional
/// \return out resulting datum
///
/// \since 1.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Count(const Datum& datum,
                    const CountOptions& options = CountOptions::Defaults(),
                    ExecContext* ctx = NULLPTR);

/// \brief Compute the mean of a numeric array.
///
/// \param[in] value datum to compute the mean, expecting Array
/// \param[in] options see ScalarAggregateOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return datum of the computed mean as a DoubleScalar
///
/// \since 1.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Mean(
    const Datum& value,
    const ScalarAggregateOptions& options = ScalarAggregateOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Compute the product of values of a numeric array.
///
/// \param[in] value datum to compute product of, expecting Array or ChunkedArray
/// \param[in] options see ScalarAggregateOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return datum of the computed sum as a Scalar
///
/// \since 6.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Product(
    const Datum& value,
    const ScalarAggregateOptions& options = ScalarAggregateOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Sum values of a numeric array.
///
/// \param[in] value datum to sum, expecting Array or ChunkedArray
/// \param[in] options see ScalarAggregateOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return datum of the computed sum as a Scalar
///
/// \since 1.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Sum(
    const Datum& value,
    const ScalarAggregateOptions& options = ScalarAggregateOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Calculate the first value of an array
///
/// \param[in] value input datum, expecting Array or ChunkedArray
/// \param[in] options see ScalarAggregateOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return datum of the computed first as Scalar
///
/// \since 13.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> First(
    const Datum& value,
    const ScalarAggregateOptions& options = ScalarAggregateOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Calculate the last value of an array
///
/// \param[in] value input datum, expecting Array or ChunkedArray
/// \param[in] options see ScalarAggregateOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return datum of the computed last as a Scalar
///
/// \since 13.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Last(
    const Datum& value,
    const ScalarAggregateOptions& options = ScalarAggregateOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Calculate the min / max of a numeric array
///
/// This function returns both the min and max as a struct scalar, with type
/// struct<min: T, max: T>, where T is the input type
///
/// \param[in] value input datum, expecting Array or ChunkedArray
/// \param[in] options see ScalarAggregateOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return resulting datum as a struct<min: T, max: T> scalar
///
/// \since 1.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> MinMax(
    const Datum& value,
    const ScalarAggregateOptions& options = ScalarAggregateOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Test whether any element in a boolean array evaluates to true.
///
/// This function returns true if any of the elements in the array evaluates
/// to true and false otherwise. Null values are ignored by default.
/// If null values are taken into account by setting ScalarAggregateOptions
/// parameter skip_nulls = false then Kleene logic is used.
/// See KleeneOr for more details on Kleene logic.
///
/// \param[in] value input datum, expecting a boolean array
/// \param[in] options see ScalarAggregateOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return resulting datum as a BooleanScalar
///
/// \since 3.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Any(
    const Datum& value,
    const ScalarAggregateOptions& options = ScalarAggregateOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Test whether all elements in a boolean array evaluate to true.
///
/// This function returns true if all of the elements in the array evaluate
/// to true and false otherwise. Null values are ignored by default.
/// If null values are taken into account by setting ScalarAggregateOptions
/// parameter skip_nulls = false then Kleene logic is used.
/// See KleeneAnd for more details on Kleene logic.
///
/// \param[in] value input datum, expecting a boolean array
/// \param[in] options see ScalarAggregateOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return resulting datum as a BooleanScalar

/// \since 3.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> All(
    const Datum& value,
    const ScalarAggregateOptions& options = ScalarAggregateOptions::Defaults(),
    ExecContext* ctx = NULLPTR);

/// \brief Calculate the modal (most common) value of a numeric array
///
/// This function returns top-n most common values and number of times they occur as
/// an array of `struct<mode: T, count: int64>`, where T is the input type.
/// Values with larger counts are returned before smaller ones.
/// If there are more than one values with same count, smaller value is returned first.
///
/// \param[in] value input datum, expecting Array or ChunkedArray
/// \param[in] options see ModeOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return resulting datum as an array of struct<mode: T, count: int64>
///
/// \since 2.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Mode(const Datum& value,
                   const ModeOptions& options = ModeOptions::Defaults(),
                   ExecContext* ctx = NULLPTR);

/// \brief Calculate the standard deviation of a numeric array
///
/// \param[in] value input datum, expecting Array or ChunkedArray
/// \param[in] options see VarianceOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return datum of the computed standard deviation as a DoubleScalar
///
/// \since 2.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Stddev(const Datum& value,
                     const VarianceOptions& options = VarianceOptions::Defaults(),
                     ExecContext* ctx = NULLPTR);

/// \brief Calculate the variance of a numeric array
///
/// \param[in] value input datum, expecting Array or ChunkedArray
/// \param[in] options see VarianceOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return datum of the computed variance as a DoubleScalar
///
/// \since 2.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Variance(const Datum& value,
                       const VarianceOptions& options = VarianceOptions::Defaults(),
                       ExecContext* ctx = NULLPTR);

/// \brief Calculate the skewness of a numeric array
///
/// \param[in] value input datum, expecting Array or ChunkedArray
/// \param[in] options see SkewOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return datum of the computed skewness as a DoubleScalar
///
/// \since 20.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Skew(const Datum& value,
                   const SkewOptions& options = SkewOptions::Defaults(),
                   ExecContext* ctx = NULLPTR);

/// \brief Calculate the kurtosis of a numeric array
///
/// \param[in] value input datum, expecting Array or ChunkedArray
/// \param[in] options see SkewOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return datum of the computed kurtosis as a DoubleScalar
///
/// \since 20.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Kurtosis(const Datum& value,
                       const SkewOptions& options = SkewOptions::Defaults(),
                       ExecContext* ctx = NULLPTR);

/// \brief Calculate the quantiles of a numeric array
///
/// \param[in] value input datum, expecting Array or ChunkedArray
/// \param[in] options see QuantileOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return resulting datum as an array
///
/// \since 4.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Quantile(const Datum& value,
                       const QuantileOptions& options = QuantileOptions::Defaults(),
                       ExecContext* ctx = NULLPTR);

/// \brief Calculate the approximate quantiles of a numeric array with T-Digest algorithm
///
/// \param[in] value input datum, expecting Array or ChunkedArray
/// \param[in] options see TDigestOptions for more information
/// \param[in] ctx the function execution context, optional
/// \return resulting datum as an array
///
/// \since 4.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> TDigest(const Datum& value,
                      const TDigestOptions& options = TDigestOptions::Defaults(),
                      ExecContext* ctx = NULLPTR);

/// \brief Find the first index of a value in an array.
///
/// \param[in] value The array to search.
/// \param[in] options The array to search for. See IndexOptions.
/// \param[in] ctx the function execution context, optional
/// \return out a Scalar containing the index (or -1 if not found).
///
/// \since 5.0.0
/// \note API not yet finalized
ARROW_EXPORT
Result<Datum> Index(const Datum& value, const IndexOptions& options,
                    ExecContext* ctx = NULLPTR);

}  // namespace compute
}  // namespace arrow
