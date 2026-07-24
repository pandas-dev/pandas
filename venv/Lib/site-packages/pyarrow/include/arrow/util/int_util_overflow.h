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
#include <initializer_list>
#include <limits>
#include <optional>
#include <type_traits>

#include "arrow/status.h"
#include "arrow/util/macros.h"
#include "arrow/util/visibility.h"

#include "arrow/vendored/safeint/safe_math.h"

namespace arrow {
namespace internal {

// Define functions AddWithOverflow, SubtractWithOverflow, MultiplyWithOverflow
// with the signature `bool(T u, T v, T* out)` where T is an integer type.
// On overflow, these functions return true.  Otherwise, false is returned
// and `out` is updated with the result of the operation.

#define SAFE_INT_OP_WITH_OVERFLOW(_func_name, _op_name, _c_type, _type)             \
  [[nodiscard]] static inline bool _func_name(_c_type u, _c_type v, _c_type* out) { \
    return !check_##_op_name##_##_type##_##_type(u, v, out);                        \
  }

#define SAFE_INT_OPS_WITH_OVERFLOW(_func_name, _op_name)            \
  SAFE_INT_OP_WITH_OVERFLOW(_func_name, _op_name, int32_t, int32)   \
  SAFE_INT_OP_WITH_OVERFLOW(_func_name, _op_name, int64_t, int64)   \
  SAFE_INT_OP_WITH_OVERFLOW(_func_name, _op_name, uint32_t, uint32) \
  SAFE_INT_OP_WITH_OVERFLOW(_func_name, _op_name, uint64_t, uint64)

SAFE_INT_OPS_WITH_OVERFLOW(SafeIntAddWithOverflow, add)
SAFE_INT_OPS_WITH_OVERFLOW(SafeIntSubtractWithOverflow, sub)
SAFE_INT_OPS_WITH_OVERFLOW(SafeIntMultiplyWithOverflow, mul)

#undef SAFE_INT_OP_WITH_OVERFLOW
#undef SAFE_INT_OPS_WITH_OVERFLOW

template <typename Int, typename SignedRet, typename UnsignedRet>
using transformed_int_t =
    std::conditional_t<std::is_signed_v<Int>, SignedRet, UnsignedRet>;

template <typename Int>
using upscaled_int32_t = transformed_int_t<Int, int32_t, uint32_t>;

// Use GCC/CLang builtins for checked arithmetic, promising better performance
// than SafeInt's hand-written implementations.
#if defined __has_builtin
#  if __has_builtin(__builtin_object_size)
#    define USE_CHECKED_ARITHMETIC_BUILTINS 1
#  else
#    define USE_CHECKED_ARITHMETIC_BUILTINS 0
#  endif
#endif

template <typename Int>
[[nodiscard]] bool AddWithOverflowGeneric(Int u, Int v, Int* out) {
#if USE_CHECKED_ARITHMETIC_BUILTINS
  return __builtin_add_overflow(u, v, out);
#else
  if constexpr (sizeof(Int) < 4) {
    using UpscaledInt = upscaled_int32_t<Int>;
    auto r = static_cast<UpscaledInt>(u) + static_cast<UpscaledInt>(v);
    *out = static_cast<Int>(r);
    return r != *out;
  } else {
    return SafeIntAddWithOverflow(u, v, out);
  }
#endif
}

template <typename Int>
[[nodiscard]] bool SubtractWithOverflowGeneric(Int u, Int v, Int* out) {
#if USE_CHECKED_ARITHMETIC_BUILTINS
  return __builtin_sub_overflow(u, v, out);
#else
  if constexpr (sizeof(Int) < 4) {
    using UpscaledInt = upscaled_int32_t<Int>;
    auto r = static_cast<UpscaledInt>(u) - static_cast<UpscaledInt>(v);
    *out = static_cast<Int>(r);
    return r != *out;
  } else {
    return SafeIntSubtractWithOverflow(u, v, out);
  }
#endif
}

template <typename Int>
[[nodiscard]] bool MultiplyWithOverflowGeneric(Int u, Int v, Int* out) {
#if USE_CHECKED_ARITHMETIC_BUILTINS
  return __builtin_mul_overflow(u, v, out);
#else
  if constexpr (sizeof(Int) < 4) {
    using UpscaledInt = upscaled_int32_t<Int>;
    auto r = static_cast<UpscaledInt>(u) * static_cast<UpscaledInt>(v);
    *out = static_cast<Int>(r);
    return r != *out;
  } else {
    return SafeIntMultiplyWithOverflow(u, v, out);
  }
#endif
}

template <typename Int>
[[nodiscard]] bool DivideWithOverflowGeneric(Int u, Int v, Int* out) {
  if (v == 0) {
    *out = Int{};
    return true;
  }
  if constexpr (std::is_signed_v<Int>) {
    constexpr auto kMin = std::numeric_limits<Int>::min();
    if (u == kMin && v == -1) {
      *out = kMin;
      return true;
    }
  }
  *out = u / v;
  return false;
}

// Define non-generic versions of the above so as to benefit from automatic
// integer conversion, to allow for mixed-type calls such as
// AddWithOverflow(int32_t, int64_t, int64_t*).

#define NON_GENERIC_OP_WITH_OVERFLOW(_func_name, _c_type)                    \
  [[nodiscard]] inline bool _func_name(_c_type u, _c_type v, _c_type* out) { \
    return ARROW_PREDICT_FALSE(_func_name##Generic(u, v, out));              \
  }

#define NON_GENERIC_OPS_WITH_OVERFLOW(_func_name)    \
  NON_GENERIC_OP_WITH_OVERFLOW(_func_name, int8_t)   \
  NON_GENERIC_OP_WITH_OVERFLOW(_func_name, uint8_t)  \
  NON_GENERIC_OP_WITH_OVERFLOW(_func_name, int16_t)  \
  NON_GENERIC_OP_WITH_OVERFLOW(_func_name, uint16_t) \
  NON_GENERIC_OP_WITH_OVERFLOW(_func_name, int32_t)  \
  NON_GENERIC_OP_WITH_OVERFLOW(_func_name, uint32_t) \
  NON_GENERIC_OP_WITH_OVERFLOW(_func_name, int64_t)  \
  NON_GENERIC_OP_WITH_OVERFLOW(_func_name, uint64_t)

NON_GENERIC_OPS_WITH_OVERFLOW(AddWithOverflow)
NON_GENERIC_OPS_WITH_OVERFLOW(SubtractWithOverflow)
NON_GENERIC_OPS_WITH_OVERFLOW(MultiplyWithOverflow)
NON_GENERIC_OPS_WITH_OVERFLOW(DivideWithOverflow)

#undef NON_GENERIC_OPS_WITH_OVERFLOW
#undef NON_GENERIC_OP_WITH_OVERFLOW

// Convenience functions over an arbitrary number of arguments
template <typename Int>
std::optional<Int> AddWithOverflow(std::initializer_list<Int> vs) {
  if (vs.size() == 0) {
    return {};
  }
  auto it = vs.begin();
  Int v = *it++;
  while (it != vs.end()) {
    if (ARROW_PREDICT_FALSE(AddWithOverflowGeneric(v, *it++, &v))) {
      return {};
    }
  }
  return v;
}

template <typename Int>
std::optional<Int> MultiplyWithOverflow(std::initializer_list<Int> vs) {
  if (vs.size() == 0) {
    return {};
  }
  auto it = vs.begin();
  Int v = *it++;
  while (it != vs.end()) {
    if (ARROW_PREDICT_FALSE(MultiplyWithOverflowGeneric(v, *it++, &v))) {
      return {};
    }
  }
  return v;
}

// Define function NegateWithOverflow with the signature `bool(T u, T* out)`
// where T is a signed integer type.  On overflow, these functions return true.
// Otherwise, false is returned and `out` is updated with the result of the
// operation.
template <typename Int>
[[nodiscard]] bool NegateWithOverflow(Int v, Int* out) {
  return SubtractWithOverflow(Int{}, v, out);
}

/// Signed addition with well-defined behaviour on overflow (as unsigned)
template <typename SignedInt>
SignedInt SafeSignedAdd(SignedInt u, SignedInt v) {
  using UnsignedInt = typename std::make_unsigned<SignedInt>::type;
  return static_cast<SignedInt>(static_cast<UnsignedInt>(u) +
                                static_cast<UnsignedInt>(v));
}

/// Signed subtraction with well-defined behaviour on overflow (as unsigned)
template <typename SignedInt>
SignedInt SafeSignedSubtract(SignedInt u, SignedInt v) {
  using UnsignedInt = typename std::make_unsigned<SignedInt>::type;
  return static_cast<SignedInt>(static_cast<UnsignedInt>(u) -
                                static_cast<UnsignedInt>(v));
}

/// Signed negation with well-defined behaviour on overflow (as unsigned)
template <typename SignedInt>
SignedInt SafeSignedNegate(SignedInt u) {
  using UnsignedInt = typename std::make_unsigned<SignedInt>::type;
  return static_cast<SignedInt>(~static_cast<UnsignedInt>(u) + 1);
}

/// Signed left shift with well-defined behaviour on negative numbers or overflow
template <typename SignedInt, typename Shift>
SignedInt SafeLeftShift(SignedInt u, Shift shift) {
  using UnsignedInt = typename std::make_unsigned<SignedInt>::type;
  return static_cast<SignedInt>(static_cast<UnsignedInt>(u) << shift);
}

}  // namespace internal
}  // namespace arrow
