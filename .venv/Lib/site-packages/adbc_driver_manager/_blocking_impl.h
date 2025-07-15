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

/// Allow KeyboardInterrupt to function with ADBC in Python.
///
/// Call SetBlockingCallback to register a callback.  This will temporarily
/// suppress the Python SIGINT handler.  When SIGINT is received, this module
/// will handle it by calling the callback.

#include <string>

namespace pyadbc_driver_manager {

/// \brief Set up internal state to handle.
/// \return An error message (or empty string).
std::string InitBlockingCallback();
/// \brief Set the callback for when SIGINT is received.
/// \return An error message (or empty string).
std::string SetBlockingCallback(void (*callback)(void*), void* data);
/// \brief Clear the callback for when SIGINT is received.
/// \return An error message (or empty string).
std::string ClearBlockingCallback();

}  // namespace pyadbc_driver_manager
