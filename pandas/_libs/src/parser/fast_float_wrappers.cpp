/*
Thin C++ wrappers exposing fast_float::from_chars with C linkage.
fast_float provides IEEE 754 round-to-even (i.e. correctly rounded) parsing.
*/
#include "fast_float/fast_float.h"

#include <system_error>

extern "C" {

int fast_float_strtod(const char *start, const char *end, double *value,
                      const char **endptr, char decimal) {
  fast_float::parse_options options{fast_float::chars_format::general, decimal};
  auto result = fast_float::from_chars_advanced(start, end, *value, options);
  // No error or overflow/underflow are valid
  if (result.ec == std::errc() || result.ec == std::errc::result_out_of_range) {
    *endptr = result.ptr;
    return 0;
  }
  return -1;
}

int try_parse_plain_double(const char *start, const char *end, char decimal,
                           double *out) {
  fast_float::parse_options options{
      fast_float::chars_format::general | fast_float::chars_format::no_infnan |
          fast_float::chars_format::allow_leading_plus,
      decimal};
  auto result = fast_float::from_chars_advanced(start, end, *out, options);
  if ((result.ec == std::errc() ||
       result.ec == std::errc::result_out_of_range) &&
      result.ptr == end) {
    return 0;
  }
  return -1;
}

} // extern "C"
