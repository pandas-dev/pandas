// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#include "pandas/util/bitarray.h"

#include <new>

#include "pandas/status.h"
#include "pandas/util.h"

namespace pandas {

BitArray::~BitArray() {
  delete[] bits_;
}

Status BitArray::Init(size_t length) {
  size_t bufsize = util::ceil_byte(length / 8);
  try {
    bits_ = new uint8_t[bufsize];
    memset(bits_, 0, bufsize);
  } catch (const std::bad_alloc& e) {
    return Status::OutOfMemory("BitArray buffer failed");
  }
  length_ = length;
  count_ = 0;
  return Status::OK();
}

} // namespace pandas
