// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#include "pandas/status.h"

#include <assert.h>
#include <stdio.h>
#include <stdint.h>

namespace pandas {


Status::Status(StatusCode code, const std::string& msg, int16_t posix_code) {
  assert(code != StatusCode::OK);
  const uint32_t size = msg.size();
  char* result = new char[size + 7];
  memcpy(result, &size, sizeof(size));
  result[4] = static_cast<char>(code);
  memcpy(result + 5, &posix_code, sizeof(posix_code));
  memcpy(result + 7, msg.c_str(), msg.size());
  state_ = result;
}

const char* Status::CopyState(const char* state) {
  uint32_t size;
  memcpy(&size, state, sizeof(size));
  char* result = new char[size + 7];
  memcpy(result, state, size + 7);
  return result;
}

std::string Status::CodeAsString() const {
  if (state_ == NULL) {
    return "OK";
  }

  const char* type;
  switch (code()) {
    case StatusCode::OK:
      type = "OK";
      break;
    case StatusCode::NotImplemented:
      type = "Not implemented";
      break;
    case StatusCode::Invalid:
      type = "Invalid argument";
      break;
    case StatusCode::KeyError:
      type = "Key error";
      break;
    case StatusCode::OutOfMemory:
      type = "Out of memory";
      break;
  }
  return std::string(type);
}

std::string Status::ToString() const {
  std::string result(CodeAsString());
  if (state_ == NULL) {
    return result;
  }

  result.append(": ");

  uint32_t length;
  memcpy(&length, state_, sizeof(length));
  result.append(reinterpret_cast<const char*>(state_ + 7), length);
  return result;
}

}  // namespace pandas
