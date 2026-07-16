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

typedef enum {
  PD_STRTOI_OK = 0,
  PD_STRTOI_OVERFLOW = 1,
  PD_STRTOI_INVALID = -1,
} pd_strtoi_status;

pd_strtoi_status pd_strtoll(const char *start, const char *end, int64_t *value,
                            const char **endptr);
pd_strtoi_status pd_strtoull(const char *start, const char *end,
                             uint64_t *value, const char **endptr);

typedef struct uint_state {
  int seen_sint;
  int seen_uint;
  int seen_null;
} uint_state;

uint64_t str_to_uint64(uint_state *state, const char *p_item, int64_t length,
                       int *error, char tsep);
int64_t str_to_int64(const char *p_item, int64_t length, int *error, char tsep);

#ifdef __cplusplus
} // extern "C"
#endif
