/*
Thin C++ wrapper exposing fast_float::from_chars with C linkage.
fast_float provides IEEE 754 round-to-even (i.e. correctly rounded) parsing.
*/
#include "fast_float/fast_float.h"

#include <system_error>

extern "C" {

int fast_float_strtod(const char *start, const char *end, double *value,
                      const char **endptr) {
  auto result = fast_float::from_chars(start, end, *value);
  if (result.ec == std::errc()) {
    *endptr = result.ptr;
    return 0;
  }
  return -1;
}

} // extern "C"
