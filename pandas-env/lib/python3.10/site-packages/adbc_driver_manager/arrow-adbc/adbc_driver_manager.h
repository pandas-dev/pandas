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

/// \file arrow-adbc/adbc_driver_manager.h ADBC Driver Manager
///
/// A helper library to dynamically load and use multiple ADBC drivers in the
/// same process.

#pragma once

#include <arrow-adbc/adbc.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef ADBC_DRIVER_MANAGER_H
#define ADBC_DRIVER_MANAGER_H

/// \brief Common entry point for drivers via the driver manager.
///
/// The driver manager can fill in default implementations of some
/// ADBC functions for drivers. Drivers must implement a minimum level
/// of functionality for this to be possible, however, and some
/// functions must be implemented by the driver.
///
/// \param[in] driver_name An identifier for the driver (e.g. a path to a
///   shared library on Linux).
/// \param[in] entrypoint An identifier for the entrypoint (e.g. the
///   symbol to call for AdbcDriverInitFunc on Linux).
/// \param[in] version The ADBC revision to attempt to initialize.
/// \param[out] driver The table of function pointers to initialize.
/// \param[out] error An optional location to return an error message
///   if necessary.
ADBC_EXPORT
AdbcStatusCode AdbcLoadDriver(const char* driver_name, const char* entrypoint,
                              int version, void* driver, struct AdbcError* error);

/// \brief Common entry point for drivers via the driver manager.
///
/// The driver manager can fill in default implementations of some
/// ADBC functions for drivers. Drivers must implement a minimum level
/// of functionality for this to be possible, however, and some
/// functions must be implemented by the driver.
///
/// \param[in] init_func The entrypoint to call.
/// \param[in] version The ADBC revision to attempt to initialize.
/// \param[out] driver The table of function pointers to initialize.
/// \param[out] error An optional location to return an error message
///   if necessary.
ADBC_EXPORT
AdbcStatusCode AdbcLoadDriverFromInitFunc(AdbcDriverInitFunc init_func, int version,
                                          void* driver, struct AdbcError* error);

/// \brief Set the AdbcDriverInitFunc to use.
///
/// This is an extension to the ADBC API. The driver manager shims
/// the AdbcDatabase* functions to allow you to specify the
/// driver/entrypoint dynamically. This function lets you set the
/// entrypoint explicitly, for applications that can dynamically
/// load drivers on their own.
ADBC_EXPORT
AdbcStatusCode AdbcDriverManagerDatabaseSetInitFunc(struct AdbcDatabase* database,
                                                    AdbcDriverInitFunc init_func,
                                                    struct AdbcError* error);

/// \brief Get a human-friendly description of a status code.
ADBC_EXPORT
const char* AdbcStatusCodeMessage(AdbcStatusCode code);

#endif  // ADBC_DRIVER_MANAGER_H

#ifdef __cplusplus
}
#endif
