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
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef ADBC_DRIVER_MANAGER_H
#define ADBC_DRIVER_MANAGER_H

typedef uint32_t AdbcLoadFlags;

#define ADBC_LOAD_FLAG_SEARCH_ENV 1
#define ADBC_LOAD_FLAG_SEARCH_USER 2
#define ADBC_LOAD_FLAG_SEARCH_SYSTEM 4
#define ADBC_LOAD_FLAG_ALLOW_RELATIVE_PATHS 8

#define ADBC_LOAD_FLAG_DEFAULT                              \
  (ADBC_LOAD_FLAG_SEARCH_ENV | ADBC_LOAD_FLAG_SEARCH_USER | \
   ADBC_LOAD_FLAG_SEARCH_SYSTEM | ADBC_LOAD_FLAG_ALLOW_RELATIVE_PATHS)

/// \brief Common entry point for drivers via the driver manager.
///
/// The driver manager can fill in default implementations of some
/// ADBC functions for drivers. Drivers must implement a minimum level
/// of functionality for this to be possible, however, and some
/// functions must be implemented by the driver.
///
/// \param[in] driver_name An identifier for the driver (e.g. a path to a
///   shared library on Linux).
/// \param[in] entrypoint An identifier for the entrypoint (e.g. the symbol to
///   call for AdbcDriverInitFunc on Linux).  If not provided, search for an
///   entrypoint based on the driver name.
/// \param[in] version The ADBC revision to attempt to initialize.
/// \param[out] driver The table of function pointers to initialize.
/// \param[out] error An optional location to return an error message
///   if necessary.
ADBC_EXPORT
AdbcStatusCode AdbcLoadDriver(const char* driver_name, const char* entrypoint,
                              int version, void* driver, struct AdbcError* error);

/// \brief Common entry point to search for and load a driver or manifest.
///
/// The driver manager can fill in default implementations of some ADBC functions
/// for drivers. Drivers must implement a minimum level of functionality for this
/// to be possible, however, and some functions must be implemented by the driver.
///
/// This function is different from AdbcLoadDriver in that it also accepts the name
/// of a driver manifest file, and allows specifying options to control what
/// directories it will search through. The behavior is as follows:
///
/// If the passed in driver_name is an absolute path:
/// - If the path has a `.toml` extension, it will attempt to parse the manifest and load
///   the driver specified within it. Erroring if this fails.
/// - If the path has an extension other than `.toml`, it will attempt to load the path as
///   a shared library. Erroring if this fails.
///
/// If the passed in driver_name does not have an extension and is not an absolute path:
/// - The load_options parameter will control whether the driver manager will search
///   the ADBC_CONFIG_PATH environment variable, the user configuration directory, and/or
///   the system level directory of /etc/adbc for either a manifest file or a shared
///   library.
/// - For each path to be searched, it will first look for <path>/<driver_name>.toml. If
///   that file exists, it will attempt to parse the manifest and load the driver
///   specified within it, erroring if this fails.
/// - If the manifest file does not exist, it will then look for
/// <path>/<driver_name>.<extension>
///   where <extension> is one of the following: `.so`, `.dll`, `.dylib`. If it can load
///   that shared library, then success is returned. Otherwise it moves to the next
///   directory until the search is either successful, or all directories have been
///   searched.
///
/// \param[in] driver_name An identifier for the driver (e.g. a path to a
///   shared library on Linux or the basename of a manifest file).
/// \param[in] entrypoint  An identifier for the entrypoint (e.g. the symbol to
///   call for AdbcDriverInitFunc on Linux).  If not provided, search for an
///   entrypoint based on the driver name.
/// \param[in] version The ADBC revision to attempt to initialize.
/// \param[in] load_options bit mask of AdbcLoadFlags to control the directories searched
/// \param[out] raw_driver The table of function pointers to initialize
/// \param[out] error An optional location to return an error message
/// @return
ADBC_EXPORT
AdbcStatusCode AdbcFindLoadDriver(const char* driver_name, const char* entrypoint,
                                  const int version, const AdbcLoadFlags load_options,
                                  void* driver, struct AdbcError* error);

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

/// \brief Set the load flags for the driver manager.
///
/// This is an extension to the ADBC API. The driver manager shims
/// the AdbcDatabase* functions to allow you to specify the
/// driver/entrypoint dynamically. This function lets you set the
/// load flags explicitly, for applications that can dynamically
/// load drivers on their own.
///
/// If this function isn't called, the default load flags are just to
/// allow relative paths, disallowing the lookups of manifests.
ADBC_EXPORT
AdbcStatusCode AdbcDriverManagerDatabaseSetLoadFlags(struct AdbcDatabase* database,
                                                     AdbcLoadFlags flags,
                                                     struct AdbcError* error);

/// \brief Get a human-friendly description of a status code.
ADBC_EXPORT
const char* AdbcStatusCodeMessage(AdbcStatusCode code);

#endif  // ADBC_DRIVER_MANAGER_H

#ifdef __cplusplus
}
#endif
