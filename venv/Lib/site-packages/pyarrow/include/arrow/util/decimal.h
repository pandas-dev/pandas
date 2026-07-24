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

#include <cstdint>
#include <iosfwd>
#include <limits>
#include <string>
#include <string_view>
#include <utility>

#include "arrow/result.h"
#include "arrow/status.h"
#include "arrow/type_fwd.h"
#include "arrow/util/basic_decimal.h"

namespace arrow {

class Decimal64;

namespace internal {

ARROW_EXPORT
Status ToArrowStatus(DecimalStatus);

}  // namespace internal

template <>
struct IntoStatus<DecimalStatus> {
  static inline Status ToStatus(DecimalStatus st) { return internal::ToArrowStatus(st); }
};

/// Represents a signed 32-bit decimal value in two's complement.
/// Calulations wrap around and overflow is ignored.
/// The max decimal precision that can be safely represented is
/// 9 significant digits.
///
/// The implementation is split into two parts :
///
/// 1. BasicDecimal32
///    - can be safely compiled to IR without references to libstdc++
/// 2. Decimal32
///    - has additional functionality on top of BasicDecimal32 to deal with
///      strings and streams
class ARROW_EXPORT Decimal32 : public BasicDecimal32 {
 public:
  /// \cond FALSE
  // (need to avoid a duplicate definition in sphinx)
  using BasicDecimal32::BasicDecimal32;
  /// \endcond

  /// \brief constructor creates a Decimal32 from a BasicDecimal32
  constexpr Decimal32(const BasicDecimal32& value) noexcept  // NOLINT runtime/explicit
      : BasicDecimal32(value) {}

  /// \brief Parse the number from a base 10 string representation
  explicit Decimal32(const std::string& value);

  /// \brief Empty constructor creates a Decimal32 with a value of 0
  /// this is required for some older compilers
  constexpr Decimal32() noexcept : BasicDecimal32() {}

  /// \brief Divide this number by right and return the result.
  ///
  /// This operation is not destructive.
  /// The answer rounds to zero. Signs work like:
  ///   21 /  5 ->  4,  1
  ///  -21 /  5 -> -4, -1
  ///   21 / -5 -> -4,  1
  ///  -21 / -5 ->  4, -1
  /// \param[in] divisor the number to divide by
  /// \return the pair of the quotient and the remainder
  Result<std::pair<Decimal32, Decimal32>> Divide(const Decimal32& divisor) const {
    std::pair<Decimal32, Decimal32> result;
    ARROW_RETURN_NOT_OK(BasicDecimal32::Divide(divisor, &result.first, &result.second));
    return result;
  }

  /// \brief Convert the Decimal32 value to a base 10 decimal string with the given scale
  std::string ToString(int32_t scale) const;

  /// \brief Convert the value to an integer string
  std::string ToIntegerString() const;

  /// \brief Cast this value to an int64_t
  explicit operator int64_t() const;

  explicit operator Decimal64() const;

  /// \brief Convert a decimal string to a Decimal value, optionally including
  /// precision and scale if they're passed in and not null.
  static Status FromString(std::string_view s, Decimal32* out, int32_t* precision,
                           int32_t* scale = NULLPTR);
  static Status FromString(const std::string& s, Decimal32* out, int32_t* precision,
                           int32_t* scale = NULLPTR);
  static Status FromString(const char* s, Decimal32* out, int32_t* precision,
                           int32_t* scale = NULLPTR);
  static Result<Decimal32> FromString(std::string_view s);
  static Result<Decimal32> FromString(const std::string& s);
  static Result<Decimal32> FromString(const char* s);

  static Result<Decimal32> FromReal(double real, int32_t precision, int32_t scale);
  static Result<Decimal32> FromReal(float real, int32_t precision, int32_t scale);

  /// \brief Convert from a big-endian byte representation. The length must be
  ///        between 1 and 4
  /// \return error status if the length is an invalid value
  static Result<Decimal32> FromBigEndian(const uint8_t* data, int32_t length);

  /// \brief Convert Decimal32 from one scale to another
  Result<Decimal32> Rescale(int32_t original_scale, int32_t new_scale) const {
    Decimal32 out;
    ARROW_RETURN_NOT_OK(BasicDecimal32::Rescale(original_scale, new_scale, &out));
    return out;
  }

  /// \brief Convert to a signed integer
  template <typename T, typename = internal::EnableIfIsOneOf<T, int32_t, int64_t>>
  Result<T> ToInteger() const {
    return static_cast<T>(value_);
  }

  /// \brief Convert to a signed integer
  template <typename T, typename = internal::EnableIfIsOneOf<T, int32_t, int64_t>>
  Status ToInteger(T* out) const {
    return ToInteger<T>().Value(out);
  }

  /// \brief Convert to a floating-point number (scaled)
  float ToFloat(int32_t scale) const;
  /// \brief Convert to a floating-point number (scaled)
  double ToDouble(int32_t scale) const;

  /// \brief Convert to a floating-point number (scaled)
  template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
  T ToReal(int32_t scale) const {
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                  "Unexpected floating-point type");
    if constexpr (std::is_same_v<T, float>) {
      return ToFloat(scale);
    } else {
      return ToDouble(scale);
    }
  }

  ARROW_FRIEND_EXPORT friend std::ostream& operator<<(std::ostream& os,
                                                      const Decimal32& decimal);
};

class ARROW_EXPORT Decimal64 : public BasicDecimal64 {
 public:
  /// \cond FALSE
  // (need to avoid a duplicate definition in sphinx)
  using BasicDecimal64::BasicDecimal64;
  /// \endcond

  /// \brief constructor creates a Decimal64 from a BasicDecimal64
  constexpr Decimal64(const BasicDecimal64& value) noexcept  // NOLINT runtime/explicit
      : BasicDecimal64(value) {}

  explicit Decimal64(const BasicDecimal32& value) noexcept
      : BasicDecimal64(static_cast<int64_t>(value.value())) {}

  /// \brief Parse the number from a base 10 string representation
  explicit Decimal64(const std::string& value);

  /// \brief Empty constructor creates a Decimal64 with a value of 0
  /// this is required for some older compilers
  constexpr Decimal64() noexcept : BasicDecimal64() {}

  /// \brief Divide this number by right and return the result.
  ///
  /// This operation is not destructive.
  /// The answer rounds to zero. Signs work like:
  ///   21 /  5 ->  4,  1
  ///  -21 /  5 -> -4, -1
  ///   21 / -5 -> -4,  1
  ///  -21 / -5 ->  4, -1
  /// \param[in] divisor the number to divide by
  /// \return the pair of the quotient and the remainder
  Result<std::pair<Decimal64, Decimal64>> Divide(const Decimal64& divisor) const {
    std::pair<Decimal64, Decimal64> result;
    ARROW_RETURN_NOT_OK(BasicDecimal64::Divide(divisor, &result.first, &result.second));
    return result;
  }

  /// \brief Convert the Decimal64 value to a base 10 decimal string with the given scale
  std::string ToString(int32_t scale) const;

  /// \brief Convert the value to an integer string
  std::string ToIntegerString() const;

  /// \brief Cast this value to an int64_t
  explicit operator int64_t() const;

  /// \brief Convert a decimal string to a Decimal value, optionally including
  /// precision and scale if they're passed in and not null.
  static Status FromString(std::string_view s, Decimal64* out, int32_t* precision,
                           int32_t* scale = NULLPTR);
  static Status FromString(const std::string& s, Decimal64* out, int32_t* precision,
                           int32_t* scale = NULLPTR);
  static Status FromString(const char* s, Decimal64* out, int32_t* precision,
                           int32_t* scale = NULLPTR);
  static Result<Decimal64> FromString(std::string_view s);
  static Result<Decimal64> FromString(const std::string& s);
  static Result<Decimal64> FromString(const char* s);

  static Result<Decimal64> FromReal(double real, int32_t precision, int32_t scale);
  static Result<Decimal64> FromReal(float real, int32_t precision, int32_t scale);

  /// \brief Convert from a big-endian byte representation. The length must be
  ///        between 1 and 4
  /// \return error status if the length is an invalid value
  static Result<Decimal64> FromBigEndian(const uint8_t* data, int32_t length);

  /// \brief Convert Decimal64 from one scale to another
  Result<Decimal64> Rescale(int32_t original_scale, int32_t new_scale) const {
    Decimal64 out;
    ARROW_RETURN_NOT_OK(BasicDecimal64::Rescale(original_scale, new_scale, &out));
    return out;
  }

  /// \brief Convert to a signed integer
  template <typename T, typename = internal::EnableIfIsOneOf<T, int32_t, int64_t>>
  Result<T> ToInteger() const {
    return static_cast<T>(value_);
  }

  /// \brief Convert to a signed integer
  template <typename T, typename = internal::EnableIfIsOneOf<T, int32_t, int64_t>>
  Status ToInteger(T* out) const {
    return ToInteger<T>().Value(out);
  }

  /// \brief Convert to a floating-point number (scaled)
  float ToFloat(int32_t scale) const;
  /// \brief Convert to a floating-point number (scaled)
  double ToDouble(int32_t scale) const;

  /// \brief Convert to a floating-point number (scaled)
  template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
  T ToReal(int32_t scale) const {
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                  "Unexpected floating-point type");
    if constexpr (std::is_same_v<T, float>) {
      return ToFloat(scale);
    } else {
      return ToDouble(scale);
    }
  }

  ARROW_FRIEND_EXPORT friend std::ostream& operator<<(std::ostream& os,
                                                      const Decimal64& decimal);
};

/// Represents a signed 128-bit integer in two's complement.
/// Calculations wrap around and overflow is ignored.
/// The max decimal precision that can be safely represented is
/// 38 significant digits.
///
/// For a discussion of the algorithms, look at Knuth's volume 2,
/// Semi-numerical Algorithms section 4.3.1.
///
/// Adapted from the Apache ORC C++ implementation
///
/// The implementation is split into two parts :
///
/// 1. BasicDecimal128
///    - can be safely compiled to IR without references to libstdc++.
/// 2. Decimal128
///    - has additional functionality on top of BasicDecimal128 to deal with
///      strings and streams.
class ARROW_EXPORT Decimal128 : public BasicDecimal128 {
 public:
  /// \cond FALSE
  // (need to avoid a duplicate definition in Sphinx)
  using BasicDecimal128::BasicDecimal128;
  /// \endcond

  /// \brief constructor creates a Decimal128 from a BasicDecimal128.
  constexpr Decimal128(const BasicDecimal128& value) noexcept  // NOLINT runtime/explicit
      : BasicDecimal128(value) {}

  /// \brief Parse the number from a base 10 string representation.
  explicit Decimal128(const std::string& value);

  /// \brief Empty constructor creates a Decimal128 with a value of 0.
  // This is required on some older compilers.
  constexpr Decimal128() noexcept : BasicDecimal128() {}

  /// Divide this number by right and return the result.
  ///
  /// This operation is not destructive.
  /// The answer rounds to zero. Signs work like:
  ///   21 /  5 ->  4,  1
  ///  -21 /  5 -> -4, -1
  ///   21 / -5 -> -4,  1
  ///  -21 / -5 ->  4, -1
  /// \param[in] divisor the number to divide by
  /// \return the pair of the quotient and the remainder
  Result<std::pair<Decimal128, Decimal128>> Divide(const Decimal128& divisor) const {
    std::pair<Decimal128, Decimal128> result;
    ARROW_RETURN_NOT_OK(BasicDecimal128::Divide(divisor, &result.first, &result.second));
    return result;
  }

  /// \brief Convert the Decimal128 value to a base 10 decimal string with the given
  /// scale.
  std::string ToString(int32_t scale) const;

  /// \brief Convert the value to an integer string
  std::string ToIntegerString() const;

  /// \brief Cast this value to an int64_t.
  explicit operator int64_t() const;

  /// \brief Convert a decimal string to a Decimal128 value, optionally including
  /// precision and scale if they're passed in and not null.
  static Status FromString(std::string_view s, Decimal128* out, int32_t* precision,
                           int32_t* scale = NULLPTR);
  static Status FromString(const std::string& s, Decimal128* out, int32_t* precision,
                           int32_t* scale = NULLPTR);
  static Status FromString(const char* s, Decimal128* out, int32_t* precision,
                           int32_t* scale = NULLPTR);
  static Result<Decimal128> FromString(std::string_view s);
  static Result<Decimal128> FromString(const std::string& s);
  static Result<Decimal128> FromString(const char* s);

  static Result<Decimal128> FromReal(double real, int32_t precision, int32_t scale);
  static Result<Decimal128> FromReal(float real, int32_t precision, int32_t scale);

  /// \brief Convert from a big-endian byte representation. The length must be
  ///        between 1 and 16.
  /// \return error status if the length is an invalid value
  static Result<Decimal128> FromBigEndian(const uint8_t* data, int32_t length);

  /// \brief Convert Decimal128 from one scale to another
  Result<Decimal128> Rescale(int32_t original_scale, int32_t new_scale) const {
    Decimal128 out;
    ARROW_RETURN_NOT_OK(BasicDecimal128::Rescale(original_scale, new_scale, &out));
    return out;
  }

  /// \brief Convert to a signed integer
  template <typename T, typename = internal::EnableIfIsOneOf<T, int32_t, int64_t>>
  Result<T> ToInteger() const {
    constexpr auto min_value = std::numeric_limits<T>::min();
    constexpr auto max_value = std::numeric_limits<T>::max();
    const auto& self = *this;
    if (self < min_value || self > max_value) {
      return Status::Invalid("Invalid cast from Decimal128 to ", sizeof(T),
                             " byte integer");
    }
    return static_cast<T>(low_bits());
  }

  /// \brief Convert to a signed integer
  template <typename T, typename = internal::EnableIfIsOneOf<T, int32_t, int64_t>>
  Status ToInteger(T* out) const {
    return ToInteger<T>().Value(out);
  }

  /// \brief Convert to a floating-point number (scaled)
  float ToFloat(int32_t scale) const;
  /// \brief Convert to a floating-point number (scaled)
  double ToDouble(int32_t scale) const;

  /// \brief Convert to a floating-point number (scaled)
  template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
  T ToReal(int32_t scale) const {
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                  "Unexpected floating-point type");
    if constexpr (std::is_same_v<T, float>) {
      return ToFloat(scale);
    } else {
      return ToDouble(scale);
    }
  }

  ARROW_FRIEND_EXPORT friend std::ostream& operator<<(std::ostream& os,
                                                      const Decimal128& decimal);
};

/// Represents a signed 256-bit integer in two's complement.
/// The max decimal precision that can be safely represented is
/// 76 significant digits.
///
/// The implementation is split into two parts :
///
/// 1. BasicDecimal256
///    - can be safely compiled to IR without references to libstdc++.
/// 2. Decimal256
///    - (TODO) has additional functionality on top of BasicDecimal256 to deal with
///      strings and streams.
class ARROW_EXPORT Decimal256 : public BasicDecimal256 {
 public:
  /// \cond FALSE
  // (need to avoid a duplicate definition in Sphinx)
  using BasicDecimal256::BasicDecimal256;
  /// \endcond

  /// \brief constructor creates a Decimal256 from a BasicDecimal256.
  constexpr Decimal256(const BasicDecimal256& value) noexcept  // NOLINT(runtime/explicit)
      : BasicDecimal256(value) {}

  /// \brief Parse the number from a base 10 string representation.
  explicit Decimal256(const std::string& value);

  /// \brief Empty constructor creates a Decimal256 with a value of 0.
  // This is required on some older compilers.
  constexpr Decimal256() noexcept : BasicDecimal256() {}

  /// \brief Convert the Decimal256 value to a base 10 decimal string with the given
  /// scale.
  std::string ToString(int32_t scale) const;

  /// \brief Convert the value to an integer string
  std::string ToIntegerString() const;

  /// \brief Convert a decimal string to a Decimal256 value, optionally including
  /// precision and scale if they're passed in and not null.
  static Status FromString(std::string_view s, Decimal256* out, int32_t* precision,
                           int32_t* scale = NULLPTR);
  static Status FromString(const std::string& s, Decimal256* out, int32_t* precision,
                           int32_t* scale = NULLPTR);
  static Status FromString(const char* s, Decimal256* out, int32_t* precision,
                           int32_t* scale = NULLPTR);
  static Result<Decimal256> FromString(std::string_view s);
  static Result<Decimal256> FromString(const std::string& s);
  static Result<Decimal256> FromString(const char* s);

  /// \brief Convert Decimal256 from one scale to another
  Result<Decimal256> Rescale(int32_t original_scale, int32_t new_scale) const {
    Decimal256 out;
    ARROW_RETURN_NOT_OK(BasicDecimal256::Rescale(original_scale, new_scale, &out));
    return out;
  }

  /// Divide this number by right and return the result.
  ///
  /// This operation is not destructive.
  /// The answer rounds to zero. Signs work like:
  ///   21 /  5 ->  4,  1
  ///  -21 /  5 -> -4, -1
  ///   21 / -5 -> -4,  1
  ///  -21 / -5 ->  4, -1
  /// \param[in] divisor the number to divide by
  /// \return the pair of the quotient and the remainder
  Result<std::pair<Decimal256, Decimal256>> Divide(const Decimal256& divisor) const {
    std::pair<Decimal256, Decimal256> result;
    ARROW_RETURN_NOT_OK(BasicDecimal256::Divide(divisor, &result.first, &result.second));
    return result;
  }

  /// \brief Convert from a big-endian byte representation. The length must be
  ///        between 1 and 32.
  /// \return error status if the length is an invalid value
  static Result<Decimal256> FromBigEndian(const uint8_t* data, int32_t length);

  static Result<Decimal256> FromReal(double real, int32_t precision, int32_t scale);
  static Result<Decimal256> FromReal(float real, int32_t precision, int32_t scale);

  /// \brief Convert to a floating-point number (scaled).
  /// May return infinity in case of overflow.
  float ToFloat(int32_t scale) const;
  /// \brief Convert to a floating-point number (scaled)
  double ToDouble(int32_t scale) const;

  /// \brief Convert to a floating-point number (scaled)
  template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
  T ToReal(int32_t scale) const {
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                  "Unexpected floating-point type");
    if constexpr (std::is_same_v<T, float>) {
      return ToFloat(scale);
    } else {
      return ToDouble(scale);
    }
  }

  ARROW_FRIEND_EXPORT friend std::ostream& operator<<(std::ostream& os,
                                                      const Decimal256& decimal);
};

/// For an integer type, return the max number of decimal digits
/// (=minimal decimal precision) it can represent.
inline Result<int32_t> MaxDecimalDigitsForInteger(Type::type type_id) {
  switch (type_id) {
    case Type::INT8:
    case Type::UINT8:
      return 3;
    case Type::INT16:
    case Type::UINT16:
      return 5;
    case Type::INT32:
    case Type::UINT32:
      return 10;
    case Type::INT64:
      return 19;
    case Type::UINT64:
      return 20;
    default:
      break;
  }
  return Status::Invalid("Not an integer type: ", type_id);
}

}  // namespace arrow
