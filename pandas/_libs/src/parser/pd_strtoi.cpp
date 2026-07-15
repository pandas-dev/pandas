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

/* Copy a numeric token without `tsep` into `output` (requires tsep != '\0').
 *
 * Returns the number of bytes written (excluding NUL), or std::nullopt if the
 * token does not fit in `output`. `*endptr` is set to the position in `str`
 * where copying stopped (the first non-digit / non-tsep char or the end of
 * input) so the caller can detect trailing garbage like "1 ," (GH#64631).
 */
static std::optional<size_t>
copy_number_without_tsep(char output[PROCESSED_WORD_CAPACITY], const char *str,
                         const char **endptr, size_t str_len, char tsep) {
  const char *p = str;
  const char *const end = str + str_len;
  char *out = output;

  if (p < end && (*p == '+' || *p == '-')) {
    *out++ = *p++;
  }

  const char *const token_end = std::find_if(
      p, end, [tsep](char ch) { return !isdigit_ascii(ch) && ch != tsep; });
  const size_t n_digits = static_cast<size_t>(token_end - p) -
                          static_cast<size_t>(std::count(p, token_end, tsep));
  if (static_cast<size_t>(out - output) + n_digits + 1 >
      PROCESSED_WORD_CAPACITY) {
    return std::nullopt;
  }

  out = std::remove_copy(p, token_end, out, tsep);
  *out = '\0';
  if (endptr != NULL) {
    *endptr = token_end;
  }
  return static_cast<size_t>(out - output);
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
    const std::optional<size_t> written =
        copy_number_without_tsep(buffer, p, &number_end, str_len, tsep);
    if (!written.has_value()) {
      // Word is too big, probably will cause an overflow
      *error = ERROR_OVERFLOW;
      return 0;
    }
    p = buffer;
    str_len = *written;
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
    const std::optional<size_t> written =
        copy_number_without_tsep(buffer, p, &number_end, str_len, tsep);
    if (!written.has_value()) {
      // Word is too big, probably will cause an overflow
      *error = ERROR_OVERFLOW;
      return 0;
    }
    p = buffer;
    str_len = *written;
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
