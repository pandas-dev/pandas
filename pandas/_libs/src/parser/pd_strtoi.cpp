/*
Thin C++ wrappers exposing std::from_chars for integer parsing with C linkage.
std::from_chars is locale-independent and skips errno, making it measurably
faster than libc strtoll/strtoull on platforms where those functions consult
the current locale on every call.
*/
#include <charconv>
#include <cstdint>
#include <system_error>

extern "C" {

// Returns:
//   0  on success (*endptr points past the last digit consumed)
//   1  on overflow / underflow
//  -1  on invalid input (no digits parsed)
int pd_strtoll(const char *start, const char *end, int64_t *value,
               const char **endptr) {
  auto result = std::from_chars(start, end, *value, 10);
  *endptr = result.ptr;
  if (result.ec == std::errc()) {
    return 0;
  }
  if (result.ec == std::errc::result_out_of_range) {
    return 1;
  }
  return -1;
}

int pd_strtoull(const char *start, const char *end, uint64_t *value,
                const char **endptr) {
  auto result = std::from_chars(start, end, *value, 10);
  *endptr = result.ptr;
  if (result.ec == std::errc()) {
    return 0;
  }
  if (result.ec == std::errc::result_out_of_range) {
    return 1;
  }
  return -1;
}

} // extern "C"
