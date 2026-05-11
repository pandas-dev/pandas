/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/parser/simd_scan.h"

#include <xsimd/xsimd.hpp>

#include <cstdint>
#include <new>

#if defined(_MSC_VER)
#  include <intrin.h>
#endif

namespace {

using batch_u8 = xsimd::batch<std::uint8_t>;
constexpr std::size_t kStep = batch_u8::size;

static_assert(kStep >= PD_SCAN_MIN_BYTES,
              "xsimd batch<uint8_t> must be at least 16 lanes wide");

static inline unsigned ctz64(std::uint64_t value) {
#if defined(_MSC_VER)
  unsigned long index;
  _BitScanForward64(&index, value);
  return static_cast<unsigned>(index);
#else
  return static_cast<unsigned>(__builtin_ctzll(value));
#endif
}

template <int N>
static inline std::size_t scan_impl(const batch_u8 *v, const char *data,
                                    std::size_t len) {
  const auto *p = reinterpret_cast<const std::uint8_t *>(data);
  std::size_t i = 0;
  for (; i + kStep <= len; i += kStep) {
    const auto chunk = batch_u8::load_unaligned(p + i);
    auto mask = (chunk == v[0]);
    for (int j = 1; j < N; ++j) {
      mask = mask | (chunk == v[j]);
    }
    if (xsimd::any(mask)) {
      return i + ctz64(mask.mask());
    }
  }
  return i;
}

} // namespace

struct pd_scanner {
  batch_u8 v[6];
  int n;
};

extern "C" {

pd_scanner *pd_scanner_create(const char *chars, int n) {
  if (n != 2 && n != 6)
    return nullptr;
  auto *scanner = new (std::nothrow) pd_scanner;
  if (!scanner)
    return nullptr;
  scanner->n = n;
  for (int j = 0; j < n; ++j) {
    scanner->v[j] = batch_u8::broadcast(static_cast<std::uint8_t>(chars[j]));
  }
  return scanner;
}

void pd_scanner_destroy(pd_scanner *scanner) { delete scanner; }

size_t pd_scanner_scan(const pd_scanner *scanner, const char *data,
                       size_t len) {
  switch (scanner->n) {
  case 2:
    return scan_impl<2>(scanner->v, data, len);
  case 6:
    return scan_impl<6>(scanner->v, data, len);
  }
  return len;
}

} // extern "C"
