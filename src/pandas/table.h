// Copyright 2015 Cloudera Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
