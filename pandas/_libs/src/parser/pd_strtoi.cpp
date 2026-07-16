/*
Integer parsing for the tokenizer, implemented in C++ (with C linkage) so the
hot paths can use fast_float::from_chars directly and have it inline.
fast_float is locale-independent, skips errno, and parses runs of 8 digits at
a time, making it measurably faster than libc strtoll/strtoull.

Behavioral contract (inherited from the tokenizer.c implementations):
- ERROR_NO_DIGITS if no digits follow the optional sign; ERROR_INVALID_CHARS
  if non-space junk follows the number; ERROR_OVERFLOW only for a clean token
  that exceeds the target type. Callers choose the uint64/float fallbacks
  based on this distinction.
- str_to_uint64 reports a leading '-' via state->seen_sint with *error = 0.
- Stripping tsep must not lose where the token ended in the original input;
  trailing junk is detected there (GH#64631).
- No heap allocation on any path.
*/
#include "pandas/parser/pd_strtoi.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <optional>
#include <span>
#include <string_view>
#include <system_error>
#include <type_traits>

#include "fast_float/fast_float.h"
#include "pandas/parser/tokenizer.h"
#include "pandas/portable.h"

// Arrow256 allows up to 76 decimal digits.
// We rounded up to the next power of 2.
constexpr size_t PROCESSED_WORD_CAPACITY = 128;

/* Copy the leading numeric token of `str` into `output`, dropping `tsep`.
 *
 * The token (optional sign, then digits/tseps) is consumed from `str`, so on
 * return `str` starts at the first unconsumed char — trailing garbage like
 * "1 ," (GH#64631) is left for the caller to detect. Returns the written
 * subspan of `output`, or std::nullopt if the token does not fit. A token
 * with no digit after the optional sign yields the sign-only subspan so the
 * subsequent from_chars reports no-digits, matching the non-tsep path
 * (",1" with tsep=',' is not 1).
 */
static std::optional<std::span<char>>
copy_number_without_tsep(std::span<char> output, std::string_view &str,
                         char tsep) {
  assert(tsep != '\0');
  auto out = output.begin();

  if (!str.empty() && (str.front() == '+' || str.front() == '-')) {
    *out++ = str.front();
    str.remove_prefix(1);
  }
  const size_t sign_len = static_cast<size_t>(out - output.begin());

  if (str.empty() || !isdigit_ascii(str.front())) {
    // No digits to copy; a leading tsep is not part of a number.
    return output.subspan(0, sign_len);
  }

  const auto token_end = std::find_if(str.begin(), str.end(), [tsep](char ch) {
    return !isdigit_ascii(ch) && ch != tsep;
  });
  const size_t token_len = static_cast<size_t>(token_end - str.begin());
  const size_t n_digits =
      token_len - static_cast<size_t>(std::count(str.begin(), token_end, tsep));
  if (sign_len + n_digits > output.size()) {
    return std::nullopt;
  }

  out = std::remove_copy(str.begin(), token_end, out, tsep);
  str.remove_prefix(token_len);
  return output.subspan(0, static_cast<size_t>(out - output.begin()));
}

/* Shared body of str_to_int64/str_to_uint64; `state` is only used for the
 * unsigned instantiation. */
template <typename T>
static T parse_integer(std::string_view str, int *error, char tsep,
                       uint_state *state = nullptr) {
  char buffer[PROCESSED_WORD_CAPACITY];
  // The token handed to from_chars: `str` itself, or its tsep-stripped copy.
  std::string_view num = str;
  // Where the number ends in `str`; on the tsep path from_chars only sees
  // the copy, so its end pointer cannot detect trailing junk (GH#64631).
  const char *number_end = nullptr;

  if (tsep != '\0' && num.find(tsep) != std::string_view::npos) {
    // copy_number_without_tsep needs the token start, so strip leading
    // spaces here; the no-tsep path leaves that to skip_white_space.
    while (!num.empty() && isspace_ascii(num.front())) {
      num.remove_prefix(1);
    }
    const std::optional<std::span<char>> copied =
        copy_number_without_tsep(buffer, num, tsep);
    if (!copied.has_value()) {
      // Word is too big, probably will cause an overflow
      *error = ERROR_OVERFLOW;
      return 0;
    }
    number_end = num.data(); // the token was consumed from `num`
    num = std::string_view(copied->data(), copied->size());
  }

  T number;
  constexpr fast_float::parse_options options(
      fast_float::chars_format::allow_leading_plus |
      fast_float::chars_format::skip_white_space);
  const auto result = fast_float::from_chars_advanced(
      num.data(), num.data() + num.size(), number, options);
  const char *const str_end = str.data() + str.size();
  const char *endptr = number_end != nullptr ? number_end : result.ptr;
  if (result.ec == std::errc::result_out_of_range) {
    // Overflow with trailing junk → INVALID_CHARS (so caller can fall through
    // to float parsing, e.g. "18446744073709551616.0"). Pure overflow →
    // OVERFLOW (caller retries as uint64).
    *error = endptr != str_end ? ERROR_INVALID_CHARS : ERROR_OVERFLOW;
    return 0;
  }
  if (result.ec != std::errc()) {
    if constexpr (!std::is_signed_v<T>) {
      // fast_float rejects '-' outright for unsigned types, leaving ptr on
      // the sign; record it for the caller's int64/uint64 reconciliation.
      if (result.ptr != num.data() + num.size() && *result.ptr == '-') {
        state->seen_sint = 1;
        *error = 0;
        return 0;
      }
    }
    // invalid_argument: no digits after the optional sign.
    *error = ERROR_NO_DIGITS;
    return 0;
  }

  // Skip trailing spaces, then require the whole input consumed.
  while (endptr != str_end && isspace_ascii(*endptr)) {
    ++endptr;
  }
  if (endptr != str_end) {
    *error = ERROR_INVALID_CHARS;
    return 0;
  }

  if constexpr (!std::is_signed_v<T>) {
    if (number > static_cast<uint64_t>(INT64_MAX)) {
      state->seen_uint = 1;
    }
  }

  *error = 0;
  return number;
}

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
