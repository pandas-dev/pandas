// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_TABLE_H
#define PANDAS_TABLE_H

#include <memory>
#include <vector>

#include "pandas/array.h"
#include "pandas/status.h"

// A table of equal-length pandas arrays
class Table {
 public:
  ArrayPtr column(size_t i) {
    return columns_[i];
  }

  size_t ncols() {
    return columns_.size();
  }

 protected:
  std::vector<ArrayPtr> columns_;
};

#endif  // PANDAS_TABLE_H
