/*
Integer parsing for the tokenizer, implemented in C++ (with C linkage) so the
hot paths can use std::from_chars and fast_float's digit helpers directly and
have them inline. std::from_chars is locale-independent and skips errno,
making it measurably faster than libc strtoll/strtoull on platforms where
those functions consult the current locale on every call.
*/
#include "pandas/parser/pd_strtoi.h"

#include <charconv>
#include <cstring>
#include <system_error>

#include "fast_float/fast_float.h"
#include "pandas/parser/tokenizer.h"
#include "pandas/portable.h"

pd_strtoi_status pd_strtoll(const char *start, const char *end, int64_t *value,
                            const char **endptr) {
  auto result = std::from_chars(start, end, *value, 10);
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
  auto result = std::from_chars(start, end, *value, 10);
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
#define PROCESSED_WORD_CAPACITY 128

/* Copy a numeric token without `tsep` into `output`.
 *
 * Returns the number of bytes written (excluding NUL), or -1 on overflow.
 * `*endptr` is set to the position in `str` where copying stopped (the first
 * non-digit / non-tsep char or the end of input) so the caller can detect
 * trailing garbage like "1 ," (GH#64631).
 */
static int copy_number_without_tsep(char output[PROCESSED_WORD_CAPACITY],
                                    const char *str, const char **endptr,
                                    size_t str_len, char tsep) {
  const char *p = str;
  const char *end = str + str_len;
  size_t bytes_written = 0;

  if (p < end && (*p == '+' || *p == '-')) {
    output[bytes_written++] = *p++;
  }

  while (p < end && (isdigit_ascii(*p) || (tsep != '\0' && *p == tsep))) {
    if (*p != tsep) {
      if (bytes_written + 1 >= PROCESSED_WORD_CAPACITY) {
        return -1;
      }
      output[bytes_written++] = *p;
    }
    p++;
  }

  output[bytes_written] = '\0';
  if (endptr != NULL) {
    *endptr = p;
  }
  return (int)bytes_written;
}

/* SWAR fast path for the common case: an optional sign followed by 1-18
 * digits filling the whole [p_item, p_item+length) span. 18 digits cannot
 * overflow int64, and a pure-digit span can contain no thousands separator,
 * spaces or trailing junk, so full consumption means the token is exactly a
 * clean integer. Digit validation and parsing use fast_float's 8-digits-at-a-
 * time helpers, which inline here (calling them from tokenizer.c instead
 * costs a per-token cross-TU call). Returns true and stores the value in
 * *result on success; anything else returns false so the caller falls
 * through to the slow path. */
static inline bool swar_try_parse_int64(const char *p_item, int64_t length,
                                        int64_t *result) {
  if (length < 1) {
    return false;
  }
  const char *p_end = p_item + length;
  const bool neg = *p_item == '-';
  const char *digits = p_item + (neg || *p_item == '+');
  size_t rem = (size_t)(p_end - digits);
  if (rem < 1 || rem > 18) {
    return false;
  }
  uint64_t acc = 0;
  while (rem >= 8) {
    const uint64_t block = fast_float::read8_to_u64(digits);
    if (!fast_float::is_made_of_eight_digits_fast(block)) {
      return false;
    }
    acc = acc * 100000000ULL + fast_float::parse_eight_digits_unrolled(block);
    digits += 8;
    rem -= 8;
  }
  for (; rem > 0; rem--, digits++) {
    const unsigned char digit = (unsigned char)(*digits - '0');
    if (digit > 9) {
      return false;
    }
    acc = acc * 10 + digit;
  }
  *result = neg ? -(int64_t)acc : (int64_t)acc;
  return true;
}

int64_t str_to_int64(const char *p_item, int64_t length, int *error,
                     char tsep) {
  int64_t swar_result;
  if (swar_try_parse_int64(p_item, length, &swar_result)) {
    *error = 0;
    return swar_result;
  }

  const char *p = p_item;
  // Skip leading spaces.
  while (isspace_ascii(*p)) {
    ++p;
  }

  // Handle sign. std::from_chars accepts '-' but rejects '+', so strip '+'
  // after verifying the following char is a digit (not another sign).
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
    const int written =
        copy_number_without_tsep(buffer, p, &number_end, str_len, tsep);
    if (written < 0) {
      // Word is too big, probably will cause an overflow
      *error = ERROR_OVERFLOW;
      return 0;
    }
    p = buffer;
    str_len = (size_t)written;
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
    const int written =
        copy_number_without_tsep(buffer, p, &number_end, str_len, tsep);
    if (written < 0) {
      // Word is too big, probably will cause an overflow
      *error = ERROR_OVERFLOW;
      return 0;
    }
    p = buffer;
    str_len = (size_t)written;
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
