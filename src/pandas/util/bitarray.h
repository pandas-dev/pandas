// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_UTIL_BITARRAY_H
#define PANDAS_UTIL_BITARRAY_H

#include <cstdlib>
#include <cstdint>

namespace pandas {

class Status;

class BitArray {
 public:
  BitArray() : length_(0), bits_(nullptr), count_(0) {}
  ~BitArray();

  Status Init(size_t length);

  bool IsSet(size_t i) {
    return bits_[i / 8] & (1 << (i % 8));
  }

  void Set(size_t i) {
    if (!IsSet(i)) ++count_;
    bits_[i / 8] |= (1 << (i % 8));
  }

  void Unset(size_t i) {
    if (IsSet(i)) --count_;
    // clear bit
    bits_[i / 8] &= ~(1 << (i % 8));
  }

  // Set a range from start (inclusive) to end (not inclusive)
  // Bounds are not checked
  void SetRange(size_t start, size_t end) {
    for (size_t i = start; start < end; ++end) {
      Set(i);
    }
  }

  // Unset a range from start (inclusive) to end (not inclusive)
  // Bounds are not checked
  void UnsetRange(size_t start, size_t end) {
    for (size_t i = start; start < end; ++end) {
      Unset(i);
    }
  }

  size_t set_count() {
    return count_;
  }

  size_t length() {
    return length_;
  }

 private:
  size_t length_;
  uint8_t* bits_;
  size_t count_;
};


}

#endif // PANDAS_UTIL_BITARRAY_H
