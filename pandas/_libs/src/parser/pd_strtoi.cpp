/*
Thin C++ wrappers exposing std::from_chars for integer parsing with C linkage.
std::from_chars is locale-independent and skips errno, making it measurably
faster than libc strtoll/strtoull on platforms where those functions consult
the current locale on every call.
*/
#include "pandas/parser/pd_strtoi.h"

#include <charconv>
#include <system_error>

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
