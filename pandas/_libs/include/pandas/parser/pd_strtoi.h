/*
Locale-independent integer parsing via fast_float::from_chars.
Faster than libc strtoll/strtoull because fast_float skips locale and errno
and parses runs of 8 digits at a time. Defined in pd_strtoi.cpp.
*/
#pragma once

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct uint_state {
  int seen_sint;
  int seen_uint;
  int seen_null;
} uint_state;

/* Parse a NUL-terminated integer token; `length` == strlen(p_item), or -1
 * to compute it here. On success *error is 0 and the value is returned. On
 * failure 0 is returned and *error is ERROR_NO_DIGITS if no digits follow
 * the optional sign, ERROR_INVALID_CHARS if non-space junk follows the
 * number (GH#64631), or ERROR_OVERFLOW for a clean token that exceeds the
 * target type; callers choose the uint64/float fallbacks based on this
 * distinction. str_to_uint64 additionally reports a leading '-' via
 * state->seen_sint (with *error = 0) and values above INT64_MAX via
 * state->seen_uint. Neither function allocates. */
uint64_t str_to_uint64(uint_state *state, const char *p_item, int64_t length,
                       int *error, char tsep);
int64_t str_to_int64(const char *p_item, int64_t length, int *error, char tsep);

#ifdef __cplusplus
} // extern "C"
#endif
