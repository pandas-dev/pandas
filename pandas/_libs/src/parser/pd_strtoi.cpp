/*
Integer parsing for the tokenizer, implemented in C++ (with C linkage) so the
hot paths can use fast_float::from_chars directly and have it inline.
fast_float is locale-independent, skips errno, and parses runs of 8 digits at
a time, making it measurably faster than libc strtoll/strtoull.
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

struct without_tsep_result {
  size_t written;  // bytes written to output, excluding the NUL
  size_t consumed; // bytes of the source consumed
};

/* Copy a numeric token without `tsep` into `output`.
 *
 * Copying stops at the first non-digit / non-tsep char or the end of input;
 * `consumed` reports how far that got so the caller can detect trailing
 * garbage like "1 ," (GH#64631). Returns std::nullopt if the token does not
 * fit in `output`.
 */
static std::optional<without_tsep_result>
copy_number_without_tsep(std::span<char> output, std::string_view str,
                         char tsep) {
  assert(tsep != '\0');
  auto out = output.begin();

  if (!str.empty() && (str.front() == '+' || str.front() == '-')) {
    *out++ = str.front();
    str.remove_prefix(1);
  }
  const size_t sign_len = static_cast<size_t>(out - output.begin());

  if (str.empty() || !isdigit_ascii(str.front())) {
    // A number must start with a digit (after any sign); in particular a
    // leading tsep is not part of a number (",1" with tsep=',' is not 1).
    // Copying no digits makes the subsequent from_chars report no-digits,
    // matching the non-tsep path.
    *out = '\0';
    return without_tsep_result{sign_len, sign_len};
  }

  const auto token_end = std::find_if(str.begin(), str.end(), [tsep](char ch) {
    return !isdigit_ascii(ch) && ch != tsep;
  });
  const size_t token_len = static_cast<size_t>(token_end - str.begin());
  const size_t n_digits =
      token_len - static_cast<size_t>(std::count(str.begin(), token_end, tsep));
  if (sign_len + n_digits + 1 > output.size()) {
    return std::nullopt;
  }

  out = std::remove_copy(str.begin(), token_end, out, tsep);
  *out = '\0';
  return without_tsep_result{static_cast<size_t>(out - output.begin()),
                             sign_len + token_len};
}

/* Shared body of str_to_int64/str_to_uint64; `state` is only used for the
 * unsigned instantiation. */
template <typename T>
static T parse_integer(const char *p_item, int64_t length, int *error,
                       char tsep, uint_state *state = nullptr) {
  const char *p = p_item;
  // length == strlen(p_item) supplied by the caller (-1 to compute here);
  // lets from_chars get its end pointer without a strlen scan.
  size_t str_len = length < 0 ? std::strlen(p) : static_cast<size_t>(length);

  char buffer[PROCESSED_WORD_CAPACITY];
  const char *number_end = nullptr;
  if (tsep != '\0' && std::memchr(p, tsep, str_len) != nullptr) {
    // copy_number_without_tsep needs the token start, so strip leading
    // spaces here; the no-tsep path leaves that to skip_white_space.
    while (isspace_ascii(*p)) {
      ++p;
      --str_len;
    }
    const std::optional<without_tsep_result> copied =
        copy_number_without_tsep(buffer, std::string_view(p, str_len), tsep);
    if (!copied.has_value()) {
      // Word is too big, probably will cause an overflow
      *error = ERROR_OVERFLOW;
      return 0;
    }
    number_end = p + copied->consumed;
    p = buffer;
    str_len = copied->written;
  }

  T number;
  constexpr fast_float::parse_options options(
      fast_float::chars_format::allow_leading_plus |
      fast_float::chars_format::skip_white_space);
  const auto result =
      fast_float::from_chars_advanced(p, p + str_len, number, options);
  const char *endptr = result.ptr;
  if (number_end != nullptr) {
    // GH#64631: detect trailing junk in the original input that
    // copy_number_without_tsep stopped at (e.g. "1 ," with tsep=',').
    endptr = number_end;
  }
  if (result.ec == std::errc::result_out_of_range) {
    // Overflow with trailing junk → INVALID_CHARS (so caller can fall through
    // to float parsing, e.g. "18446744073709551616.0"). Pure overflow (endptr
    // at NUL) → OVERFLOW (caller retries as uint64).
    *error = *endptr ? ERROR_INVALID_CHARS : ERROR_OVERFLOW;
    return 0;
  }
  if (result.ec != std::errc()) {
    if constexpr (!std::is_signed_v<T>) {
      // fast_float rejects '-' outright for unsigned types (leaving ptr on
      // the sign); record it for the caller's int64/uint64 reconciliation.
      if (*result.ptr == '-') {
        state->seen_sint = 1;
        *error = 0;
        return 0;
      }
    }
    // invalid_argument: no digits after the optional sign.
    *error = ERROR_NO_DIGITS;
    return 0;
  }

  // Skip trailing spaces.
  while (isspace_ascii(*endptr)) {
    ++endptr;
  }

  // Did we use up all the characters?
  if (*endptr) {
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
  return parse_integer<int64_t>(p_item, length, error, tsep);
}

uint64_t str_to_uint64(uint_state *state, const char *p_item, int64_t length,
                       int *error, char tsep) {
  return parse_integer<uint64_t>(p_item, length, error, tsep, state);
}
