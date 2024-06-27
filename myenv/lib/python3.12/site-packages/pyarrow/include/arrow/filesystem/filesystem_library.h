// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#pragma once

#include "arrow/filesystem/filesystem.h"

namespace arrow::fs {
extern "C" {

// ARROW_FORCE_EXPORT ensures this function's visibility is
// _declspec(dllexport)/[[gnu::visibility("default")]] even when
// this header is #included by a non-arrow source, as in a third
// party filesystem implementation.
ARROW_FORCE_EXPORT void* arrow_filesystem_get_registry() {
  // In the case where libarrow is linked statically both to the executable and to a
  // dynamically loaded filesystem implementation library, the library contains a
  // duplicate definition of the registry into which the library's instances of
  // FileSystemRegistrar insert their factories. This function is made accessible to
  // dlsym/GetProcAddress to enable detection of such duplicate registries and merging
  // into the registry accessible to the executable.
  return internal::GetFileSystemRegistry();
}
}
}  // namespace arrow::fs
