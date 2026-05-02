/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/simd/moments_simd.hpp"

namespace pandas::moments {

template Moments accumulate_moments_simd::operator()<xsimd::avx512cd>(
    xsimd::avx512cd, std::span<const double>, bool,
    std::optional<std::span<const uint8_t>>, int) noexcept;

} // namespace pandas::moments
