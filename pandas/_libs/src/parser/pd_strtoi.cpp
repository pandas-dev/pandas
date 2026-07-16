/*
Integer parsing for the tokenizer, implemented in C++ (with C linkage) so the
hot paths can use fast_float::from_chars directly and have it inline.
fast_float::from_chars is locale-independent, skips errno, and parses runs of
8 digits at a time, making it measurably faster than libc strtoll/strtoull.
*/
#include "pandas/parser/pd_strtoi.h"

#include <algorithm>
#include <cstring>
#include <optional>
#include <span>
#include <string_view>
#include <system_error>

#include "fast_float/fast_float.h"
#include "pandas/parser/tokenizer.h"
#include "pandas/portable.h"

pd_strtoi_status pd_strtoll(const char *start, const char *end, int64_t *value,
                            const char **endptr) {
  auto result = fast_float::from_chars(start, end, *value, 10);
  *endptr = result.ptr;
  if (result.ec == std::errc()) {
    return PD_STRTOI_OK;
  }
  if (result.ec == std::errc::result_out_of_range) {
    return PD_STRTOI_OVERFLOW;
  }
  return PD_STRTOI_INVALID;
}

pd_strtoi_status pd_strtoull(const char *start, const char *end,
                             uint64_t *value, const char **endptr) {
  auto result = fast_float::from_chars(start, end, *value, 10);
  *endptr = result.ptr;
  if (result.ec == std::errc()) {
    return PD_STRTOI_OK;
  }
  if (result.ec == std::errc::result_out_of_range) {
    return PD_STRTOI_OVERFLOW;
  }
  return PD_STRTOI_INVALID;
}

// Arrow256 allows up to 76 decimal digits.
// We rounded up to the next power of 2.
constexpr size_t PROCESSED_WORD_CAPACITY = 128;

struct without_tsep_result {
  size_t written;  // bytes written to output, excluding the NUL
  size_t consumed; // bytes of the source consumed
};

/* Copy a numeric token without `tsep` into `output` (requires tsep != '\0').
 *
 * Copying stops at the first non-digit / non-tsep char or the end of input;
 * `consumed` reports how far that got so the caller can detect trailing
 * garbage like "1 ," (GH#64631). Returns std::nullopt if the token does not
 * fit in `output`.
 */
static std::optional<without_tsep_result>
copy_number_without_tsep(std::span<char> output, std::string_view str,
                         char tsep) {
  auto out = output.begin();
  std::string_view rest = str;

  if (!rest.empty() && (rest.front() == '+' || rest.front() == '-')) {
    *out++ = rest.front();
    rest.remove_prefix(1);
  }

  const auto token_end =
      std::find_if(rest.begin(), rest.end(), [tsep](char ch) {
        return !isdigit_ascii(ch) && ch != tsep;
      });
  const size_t token_len = static_cast<size_t>(token_end - rest.begin());
  const size_t n_digits =
      token_len -
      static_cast<size_t>(std::count(rest.begin(), token_end, tsep));
  const size_t sign_len = static_cast<size_t>(out - output.begin());
  if (sign_len + n_digits + 1 > output.size()) {
    return std::nullopt;
  }

  out = std::remove_copy(rest.begin(), token_end, out, tsep);
  *out = '\0';
  return without_tsep_result{static_cast<size_t>(out - output.begin()),
                             sign_len + token_len};
}

int64_t str_to_int64(const char *p_item, int64_t length, int *error,
                     char tsep) {
  const char *p = p_item;
  // Skip leading spaces.
  while (isspace_ascii(*p)) {
    ++p;
  }

  // Handle sign. fast_float::from_chars accepts '-' but rejects '+', so strip
  // '+' after verifying the following char is a digit (not another sign).
  const bool has_sign = *p == '-' || *p == '+';
  const char *digit_start = has_sign ? p + 1 : p;
  if (!isdigit_ascii(*digit_start)) {
    *error = ERROR_NO_DIGITS;
    return 0;
  }
  if (*p == '+') {
    ++p;
  }

  char buffer[PROCESSED_WORD_CAPACITY];
  // length == strlen(p_item) supplied by the caller (-1 to compute here);
  // lets from_chars get its end pointer without a strlen scan.
  size_t str_len =
      length < 0 ? strlen(p) : (size_t)length - (size_t)(p - p_item);
  const char *number_end = NULL;
  if (tsep != '\0' && memchr(p, tsep, str_len) != NULL) {
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

  int64_t number;
  const char *endptr;
  const pd_strtoi_status status = pd_strtoll(p, p + str_len, &number, &endptr);
  if (number_end != NULL) {
    // GH#64631: detect trailing junk in the original input that
    // copy_number_without_tsep stopped at (e.g. "1 ," with tsep=',').
    endptr = number_end;
  }
  if (status == PD_STRTOI_OVERFLOW) {
    // Overflow with trailing junk → INVALID_CHARS (so caller can fall through
    // to float parsing, e.g. "18446744073709551616.0"). Pure overflow (endptr
    // at NUL) → OVERFLOW (caller retries as uint64).
    *error = *endptr ? ERROR_INVALID_CHARS : ERROR_OVERFLOW;
    return 0;
  }
  if (status == PD_STRTOI_INVALID) {
    *error = ERROR_INVALID_CHARS;
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

  *error = 0;
  return number;
}

uint64_t str_to_uint64(uint_state *state, const char *p_item, int64_t length,
                       int *error, char tsep) {
  const char *p = p_item;
  // Skip leading spaces.
  while (isspace_ascii(*p)) {
    ++p;
  }

  // Handle sign.
  if (*p == '-') {
    state->seen_sint = 1;
    *error = 0;
    return 0;
  } else if (*p == '+') {
    p++;
  }

  // Check that there is a first digit.
  if (!isdigit_ascii(*p)) {
    *error = ERROR_NO_DIGITS;
    return 0;
  }

  char buffer[PROCESSED_WORD_CAPACITY];
  // length == strlen(p_item) supplied by the caller (-1 to compute here);
  // lets from_chars get its end pointer without a strlen scan.
  size_t str_len =
      length < 0 ? strlen(p) : (size_t)length - (size_t)(p - p_item);
  const char *number_end = NULL;
  if (tsep != '\0' && memchr(p, tsep, str_len) != NULL) {
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

  uint64_t number;
  const char *endptr;
  const pd_strtoi_status status = pd_strtoull(p, p + str_len, &number, &endptr);
  if (number_end != NULL) {
    // GH#64631: detect trailing junk in the original input.
    endptr = number_end;
  }
  if (status == PD_STRTOI_OVERFLOW) {
    *error = *endptr ? ERROR_INVALID_CHARS : ERROR_OVERFLOW;
    return 0;
  }
  if (status == PD_STRTOI_INVALID) {
    *error = ERROR_INVALID_CHARS;
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

  if (number > (uint64_t)INT64_MAX) {
    state->seen_uint = 1;
  }

  *error = 0;
  return number;
}
