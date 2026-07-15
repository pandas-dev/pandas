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

#ifdef __cplusplus
} // extern "C"
#endif
