/*
Integer parsing for the tokenizer, implemented in C++ (with C linkage) so the
hot paths can use fast_float::from_chars directly and have it inline.
fast_float is locale-independent, skips errno, and parses runs of 8 digits at
a time, making it measurably faster than libc strtoll/strtoull.
See pd_strtoi.h for the error-code contract these functions preserve.
*/
#include "pandas/parser/pd_strtoi.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <cstring>
#include <span>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>

#include "fast_float/fast_float.h"
#include "pandas/parser/tokenizer.h"

namespace {

// Arrow256 allows up to 76 decimal digits.
// We rounded up to the next power of 2.
constexpr size_t processed_word_capacity = 128;

// The ASCII whitespace accepted around numbers; matches both isspace_ascii
// and fast_float::is_space (relied on by skip_white_space below).
constexpr std::string_view ascii_spaces = " \t\n\v\f\r";

// Result of copy_number_without_tsep, following the from_chars conventions:
// `token` views the stripped copy in the caller's buffer, `ptr` is the first
// unconsumed char of the source, and `ec` is errc::value_too_large if the
// token does not fit in the buffer.
struct without_tsep_result {
  std::string_view token;
  const char *ptr;
  std::errc ec;
};

/* Copy the leading numeric token of `str` (optional sign, then digits/tseps)
 * into `output`, dropping `tsep`.
 *
 * Trailing garbage like "1 ," (GH#64631) is left past `ptr` for the caller
 * to detect. A token with no digit after the optional sign yields a
 * sign-only `token` so the subsequent from_chars reports no-digits,
 * matching the non-tsep path (",1" with tsep=',' is not 1).
 */
without_tsep_result copy_number_without_tsep(std::span<char> output,
                                             std::string_view str, char tsep) {
  assert(tsep != '\0');
  auto out = output.begin();

  std::ptrdiff_t sign_len = 0;
  if (!str.empty() && (str.front() == '+' || str.front() == '-')) {
    *out++ = str.front();
    str.remove_prefix(1);
    sign_len = 1;
  }

  if (str.empty() || !std::isdigit(static_cast<unsigned char>(str.front()))) {
    // No digits to copy; a leading tsep is not part of a number.
    return {.token = std::string_view(output.data(), sign_len),
            .ptr = str.data()};
  }

  auto it = str.begin();
  for (; it != str.end() &&
         (std::isdigit(static_cast<unsigned char>(*it)) || *it == tsep);
       ++it) {
    if (*it == tsep) {
      continue;
    }
    if (out == output.end()) {
      return {.ptr = str.data(), .ec = std::errc::value_too_large};
    }
    *out++ = *it;
  }
  return {.token = std::string_view(output.data(), out - output.begin()),
          .ptr = str.data() + (it - str.begin())};
}

/* Shared body of str_to_int64/str_to_uint64. */
template <typename T>
T parse_integer(std::string_view str, int *error, char tsep,
                uint_state *state = nullptr) {
  if constexpr (std::is_unsigned_v<T>) {
    assert(state != nullptr); // for seen_sint/seen_uint
  }
  // `str` is consumed while parsing; boundary checks use the original end.
  const char *const str_end = str.data() + str.size();
  std::array<char, processed_word_capacity> buffer{};
  // Number end in the original input; the tsep path parses a copy, so
  // from_chars' own end pointer cannot flag trailing junk (GH#64631).
  const char *number_end = nullptr;

  if (tsep != '\0' && str.find(tsep) != std::string_view::npos) {
    // copy_number_without_tsep needs the token start, so strip leading
    // spaces here; the no-tsep path leaves that to skip_white_space.
    str.remove_prefix(
        std::min(str.find_first_not_of(ascii_spaces), str.size()));
    const without_tsep_result copied =
        copy_number_without_tsep(buffer, str, tsep);
    if (copied.ec != std::errc()) {
      // Word is too big, probably will cause an overflow
      *error = ERROR_OVERFLOW;
      return 0;
    }
    number_end = copied.ptr;
    str = copied.token;
  }

  T number;
  constexpr fast_float::parse_options options(
      fast_float::chars_format::allow_leading_plus |
      fast_float::chars_format::skip_white_space);
  const auto result = fast_float::from_chars_advanced(
      str.data(), str.data() + str.size(), number, options);
  const char *endptr = number_end != nullptr ? number_end : result.ptr;
  if (result.ec == std::errc::result_out_of_range) {
    // Overflow with trailing junk → INVALID_CHARS (so caller can fall through
    // to float parsing, e.g. "18446744073709551616.0"). Pure overflow →
    // OVERFLOW (caller retries as uint64).
    *error = endptr != str_end ? ERROR_INVALID_CHARS : ERROR_OVERFLOW;
    return 0;
  }
  if (result.ec != std::errc()) {
    if constexpr (std::is_unsigned_v<T>) {
      // fast_float rejects '-' outright for unsigned types, leaving ptr on
      // the sign; record it for the caller's int64/uint64 reconciliation.
      // On an empty token ptr == end, and dereferencing it would read an
      // uninitialized buffer byte — hence the first check.
      if (result.ptr != str.data() + str.size() && *result.ptr == '-') {
        state->seen_sint = 1;
        *error = 0;
        return 0;
      }
    }
    // invalid_argument: no digits after the optional sign.
    *error = ERROR_NO_DIGITS;
    return 0;
  }

  // Only trailing spaces may follow the number.
  const std::string_view trailing(endptr, str_end - endptr);
  if (trailing.find_first_not_of(ascii_spaces) != std::string_view::npos) {
    *error = ERROR_INVALID_CHARS;
    return 0;
  }

  if constexpr (std::is_unsigned_v<T>) {
    if (std::cmp_greater(number, INT64_MAX)) {
      state->seen_uint = 1;
    }
  }

  *error = 0;
  return number;
}

} // namespace

int64_t str_to_int64(const char *p_item, int64_t length, int *error,
                     char tsep) {
  // length == strlen(p_item) supplied by the caller (-1 to compute here);
  // lets from_chars get its end pointer without a strlen scan.
  const size_t str_len =
      length < 0 ? std::strlen(p_item) : static_cast<size_t>(length);
  return parse_integer<int64_t>(std::string_view(p_item, str_len), error, tsep);
}

uint64_t str_to_uint64(uint_state *state, const char *p_item, int64_t length,
                       int *error, char tsep) {
  const size_t str_len =
      length < 0 ? std::strlen(p_item) : static_cast<size_t>(length);
  return parse_integer<uint64_t>(std::string_view(p_item, str_len), error, tsep,
                                 state);
}
