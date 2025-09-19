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

#include <atomic>
#include <cstdint>
#include <optional>
#include <thread>
#include <unordered_map>
#include <vector>

#include "arrow/compute/expression.h"
#include "arrow/compute/type_fwd.h"
#include "arrow/result.h"
#include "arrow/util/cpu_info.h"
#include "arrow/util/simd.h"

#if defined(__clang__) || defined(__GNUC__)
#  define BYTESWAP(x) __builtin_bswap64(x)
#  define ROTL(x, n) (((x) << (n)) | ((x) >> ((-n) & 31)))
#  define ROTL64(x, n) (((x) << (n)) | ((x) >> ((-n) & 63)))
#elif defined(_MSC_VER)
#  include <intrin.h>
#  define BYTESWAP(x) _byteswap_uint64(x)
#  define ROTL(x, n) _rotl((x), (n))
#  define ROTL64(x, n) _rotl64((x), (n))
#endif

namespace arrow {
namespace util {

// Some platforms typedef int64_t as long int instead of long long int,
// which breaks the _mm256_i64gather_epi64 and _mm256_i32gather_epi64 intrinsics
// which need long long.
// We use the cast to the type below in these intrinsics to make the code
// compile in all cases.
//
using int64_for_gather_t = const long long int;  // NOLINT runtime-int

// All MiniBatch... classes use TempVectorStack for vector allocations and can
// only work with vectors up to 1024 elements.
//
// They should only be allocated on the stack to guarantee the right sequence
// of allocation and deallocation of vectors from TempVectorStack.
//
class MiniBatch {
 public:
  static constexpr int kLogMiniBatchLength = 10;
  static constexpr int kMiniBatchLength = 1 << kLogMiniBatchLength;
};

namespace bit_util {

ARROW_EXPORT void bits_to_indexes(int bit_to_search, int64_t hardware_flags,
                                  const int num_bits, const uint8_t* bits,
                                  int* num_indexes, uint16_t* indexes,
                                  int bit_offset = 0);

ARROW_EXPORT void bits_filter_indexes(int bit_to_search, int64_t hardware_flags,
                                      const int num_bits, const uint8_t* bits,
                                      const uint16_t* input_indexes, int* num_indexes,
                                      uint16_t* indexes, int bit_offset = 0);

// Input and output indexes may be pointing to the same data (in-place filtering).
ARROW_EXPORT void bits_split_indexes(int64_t hardware_flags, const int num_bits,
                                     const uint8_t* bits, int* num_indexes_bit0,
                                     uint16_t* indexes_bit0, uint16_t* indexes_bit1,
                                     int bit_offset = 0);

// Bit 1 is replaced with byte 0xFF.
ARROW_EXPORT void bits_to_bytes(int64_t hardware_flags, const int num_bits,
                                const uint8_t* bits, uint8_t* bytes, int bit_offset = 0);

// Return highest bit of each byte.
ARROW_EXPORT void bytes_to_bits(int64_t hardware_flags, const int num_bits,
                                const uint8_t* bytes, uint8_t* bits, int bit_offset = 0);

ARROW_EXPORT bool are_all_bytes_zero(int64_t hardware_flags, const uint8_t* bytes,
                                     uint32_t num_bytes);

#if defined(ARROW_HAVE_RUNTIME_AVX2) && defined(ARROW_HAVE_RUNTIME_BMI2)
// The functions below use BMI2 instructions, be careful before calling!

namespace avx2 {
ARROW_EXPORT void bits_filter_indexes_avx2(int bit_to_search, const int num_bits,
                                           const uint8_t* bits,
                                           const uint16_t* input_indexes,
                                           int* num_indexes, uint16_t* indexes);
ARROW_EXPORT void bits_to_indexes_avx2(int bit_to_search, const int num_bits,
                                       const uint8_t* bits, int* num_indexes,
                                       uint16_t* indexes, uint16_t base_index = 0);
ARROW_EXPORT void bits_to_bytes_avx2(const int num_bits, const uint8_t* bits,
                                     uint8_t* bytes);
ARROW_EXPORT void bytes_to_bits_avx2(const int num_bits, const uint8_t* bytes,
                                     uint8_t* bits);
ARROW_EXPORT bool are_all_bytes_zero_avx2(const uint8_t* bytes, uint32_t num_bytes);
}  // namespace avx2

#endif

}  // namespace bit_util
}  // namespace util

namespace compute {

/// Modify an Expression with pre-order and post-order visitation.
/// `pre` will be invoked on each Expression. `pre` will visit Calls before their
/// arguments, `post_call` will visit Calls (and no other Expressions) after their
/// arguments. Visitors should return the Identical expression to indicate no change; this
/// will prevent unnecessary construction in the common case where a modification is not
/// possible/necessary/...
///
/// If an argument was modified, `post_call` visits a reconstructed Call with the modified
/// arguments but also receives a pointer to the unmodified Expression as a second
/// argument. If no arguments were modified the unmodified Expression* will be nullptr.
template <typename PreVisit, typename PostVisitCall>
Result<Expression> ModifyExpression(Expression expr, const PreVisit& pre,
                                    const PostVisitCall& post_call) {
  ARROW_ASSIGN_OR_RAISE(expr, Result<Expression>(pre(std::move(expr))));

  auto call = expr.call();
  if (!call) return expr;

  bool at_least_one_modified = false;
  std::vector<Expression> modified_arguments;

  for (size_t i = 0; i < call->arguments.size(); ++i) {
    ARROW_ASSIGN_OR_RAISE(auto modified_argument,
                          ModifyExpression(call->arguments[i], pre, post_call));

    if (Identical(modified_argument, call->arguments[i])) {
      continue;
    }

    if (!at_least_one_modified) {
      modified_arguments = call->arguments;
      at_least_one_modified = true;
    }

    modified_arguments[i] = std::move(modified_argument);
  }

  if (at_least_one_modified) {
    // reconstruct the call expression with the modified arguments
    auto modified_call = *call;
    modified_call.arguments = std::move(modified_arguments);
    return post_call(Expression(std::move(modified_call)), &expr);
  }

  return post_call(std::move(expr), NULLPTR);
}

// Helper class to calculate the modified number of rows to process using SIMD.
//
// Some array elements at the end will be skipped in order to avoid buffer
// overrun, when doing memory loads and stores using larger word size than a
// single array element.
//
class TailSkipForSIMD {
 public:
  static int64_t FixBitAccess(int num_bytes_accessed_together, int64_t num_rows,
                              int bit_offset) {
    int64_t num_bytes = bit_util::BytesForBits(num_rows + bit_offset);
    int64_t num_bytes_safe =
        std::max(static_cast<int64_t>(0LL), num_bytes - num_bytes_accessed_together + 1);
    int64_t num_rows_safe =
        std::max(static_cast<int64_t>(0LL), 8 * num_bytes_safe - bit_offset);
    return std::min(num_rows_safe, num_rows);
  }
  static int64_t FixBinaryAccess(int num_bytes_accessed_together, int64_t num_rows,
                                 int64_t length) {
    int64_t num_rows_to_skip = bit_util::CeilDiv(length, num_bytes_accessed_together);
    int64_t num_rows_safe =
        std::max(static_cast<int64_t>(0LL), num_rows - num_rows_to_skip);
    return num_rows_safe;
  }
  static int64_t FixVarBinaryAccess(int num_bytes_accessed_together, int64_t num_rows,
                                    const uint32_t* offsets) {
    // Do not process rows that could read past the end of the buffer using N
    // byte loads/stores.
    //
    int64_t num_rows_safe = num_rows;
    while (num_rows_safe > 0 &&
           offsets[num_rows_safe] + num_bytes_accessed_together > offsets[num_rows]) {
      --num_rows_safe;
    }
    return num_rows_safe;
  }
  static int FixSelection(int64_t num_rows_safe, int num_selected,
                          const uint16_t* selection) {
    int num_selected_safe = num_selected;
    while (num_selected_safe > 0 && selection[num_selected_safe - 1] >= num_rows_safe) {
      --num_selected_safe;
    }
    return num_selected_safe;
  }
};

}  // namespace compute
}  // namespace arrow
