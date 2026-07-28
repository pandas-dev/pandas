/*
Thin C++ wrappers exposing fast_float::from_chars for integer parsing with C
linkage. Like std::from_chars it is locale-independent and skips errno, and it
follows the same (value, ptr, ec) contract, but it consumes runs of eight
digits at a time instead of one digit per iteration.
*/
#include "pandas/parser/pd_strtoi.h"

#include <system_error>

#include "fast_float/fast_float.h"

namespace {

// fast_float::from_chars and std::from_chars agree on the surface used here:
// no leading '+', no leading whitespace, '-' rejected for unsigned types, and
// on overflow `ptr` still points past the final digit so the caller can tell a
// pure overflow from one with trailing junk.
template <typename T>
pd_strtoi_status from_chars_to_status(const char *start, const char *end,
                                      T *value, const char **endptr) {
  const auto result = fast_float::from_chars(start, end, *value, 10);
  *endptr = result.ptr;
  if (result.ec == std::errc()) {
    return PD_STRTOI_OK;
  }
  if (result.ec == std::errc::result_out_of_range) {
    return PD_STRTOI_OVERFLOW;
  }
  return PD_STRTOI_INVALID;
}

} // namespace

pd_strtoi_status pd_strtoll(const char *start, const char *end, int64_t *value,
                            const char **endptr) {
  return from_chars_to_status(start, end, value, endptr);
}

pd_strtoi_status pd_strtoull(const char *start, const char *end,
                             uint64_t *value, const char **endptr) {
  return from_chars_to_status(start, end, value, endptr);
}
