/*
Thin C++ wrapper exposing fast_float::from_chars with C linkage.
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

/* Hot-path double parse for read_csv with default settings (no thousands
 * separator, 'e'/'E' exponent). Succeeds (returns 0, writes *out) only
 * when a plain numeric token parses cleanly and consumes exactly
 * [start, end); anything else — leading/trailing spaces, leading '+',
 * inf/nan spellings, junk — returns nonzero so the caller retries through
 * the full converter with its legacy semantics.
 */
int pd_fast_double(const char *start, const char *end, char decimal,
                   double *out) {
  fast_float::parse_options options{fast_float::chars_format::general |
                                        fast_float::chars_format::no_infnan,
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
