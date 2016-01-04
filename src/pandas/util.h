// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_UTIL_H
#define PANDAS_UTIL_H

#include <cstdlib>
#include <cstring>

#define DISALLOW_COPY_AND_ASSIGN(TypeName)      \
  TypeName(const TypeName&) = delete;           \
  void operator=(const TypeName&) = delete

namespace pandas {

namespace util {

static inline size_t ceil_byte(size_t size) {
  return (size + 7) & ~7;
}

static inline size_t ceil_2bytes(size_t size) {
  return (size + 15) & ~15;
}

static inline bool get_bit(const uint8_t* bits, size_t i) {
  return bits[i / 8] & (1 << (i % 8));
}

static inline void set_bit(uint8_t* bits, size_t i, bool is_set) {
  bits[i / 8] |= (1 << (i % 8)) * is_set;
}


static inline size_t next_power2(size_t n) {
  n--;
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  if (sizeof(size_t) == 8) {
    n |= n >> 32;
  }
  n++;
  return n;
}

static void bytes_to_bits(uint8_t* bytes, size_t length, uint8_t* bits) {
  for (size_t i = 0; i < length; ++i) {
    set_bit(bits, i, static_cast<bool>(bytes[i]));
  }
}

static uint8_t* bytes_to_bits(uint8_t* bytes, size_t length, size_t* out_length) {
  if (!length) {
    return nullptr;
  }
  size_t bit_length = *out_length = ceil_byte(length) / 8;

  // TODO: it would be better to do this in preallocated memory
  uint8_t* result = reinterpret_cast<uint8_t*>(malloc(bit_length));

  if (result == nullptr) {
    // malloc failed
    return result;
  }

  memset(result, 0, bit_length);
  bytes_to_bits(bytes, length, result);
  return result;
}

} // namespace util

} // namespace pandas

#endif // PANDAS_UTIL_H
