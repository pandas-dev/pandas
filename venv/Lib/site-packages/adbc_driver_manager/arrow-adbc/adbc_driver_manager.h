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
/// - The load_options parameter will control whether the driver manager will search the
///   environment variable ADBC_DRIVER_PATH and (if built or installed with conda) the
///   conda environment, the user-level configuration, and/or the system-level
///   configuration for either a manifest file or a shared library.
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
/// \param[in] additional_search_path_list A list of additional paths to search for
///    delimited by the OS specific path list separator.
/// \param[out] driver The table of function pointers to initialize
/// \param[out] error An optional location to return an error message
ADBC_EXPORT
AdbcStatusCode AdbcFindLoadDriver(const char* driver_name, const char* entrypoint,
                                  const int version, const AdbcLoadFlags load_options,
                                  const char* additional_search_path_list, void* driver,
                                  struct AdbcError* error);

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

/// \brief Set an additional manifest search path list for the driver manager.
///
/// This is an extension to the ADBC API. The driver manager shims
/// the AdbcDatabase* functions to allow you to specify the
/// driver/entrypoint dynamically. This function lets you explicitly
/// set a path list at runtime for additional paths to search when
/// looking for a driver manifest. While users can add additional
/// paths via the ADBC_DRIVER_PATH environment variable, this allows
/// an application to specify search paths at runtime which are not tied
/// to the load flags.
///
/// Calling this function with NULL as the `path_list` will clear any
/// previously set additional search paths.
ADBC_EXPORT
AdbcStatusCode AdbcDriverManagerDatabaseSetAdditionalSearchPathList(
    struct AdbcDatabase* database, const char* path_list, struct AdbcError* error);

/// \brief Get a human-friendly description of a status code.
ADBC_EXPORT
const char* AdbcStatusCodeMessage(AdbcStatusCode code);

/// \defgroup adbc-driver-manager-connection-profile Connection Profiles
/// The ADBC driver manager can support "connection profiles" that specify
/// a driver and options to use when connecting. This allows users to
/// specify connection information in a file or environment variable, and have the
/// driver manager load the appropriate driver and set options accordingly.
///
/// This allows creating reusable connection configurations for sharing and distribution
/// without needing to hardcode driver names and options in application code. Profiles
/// will be loaded during DatabaseInit before attempting to initialize the driver. Any
/// options specified by the profile will be applied but will not override options
/// that have already been set using DatabaseSetOption.
///
/// To facilitate customization, we define an interface for implementing a Connection
/// Profile object along with a provider function definition which can be set into
/// the driver manager to allow for customized profile loading.
///
/// A profile can be specified to the Driver Manager in one of two ways,
/// which will invoke the profile provider during the call to DatabaseInit:
///
/// 1. The "profile" option can be set using DatabaseSetOption with the name of the
/// profile to load.
/// 2. The "uri" being used can have the form "profile://<profile>"
///
/// @{

/// \brief Abstract interface for connection profile providers
struct ADBC_EXPORT AdbcConnectionProfile {
  /// \brief Opaque implementation-defined state.
  /// This field is NULL if the profile is uninitialized/freed (but
  /// it need not have a value even if the profile is initialized).
  void* private_data;

  /// \brief Release the profile and perform any cleanup.
  void (*release)(struct AdbcConnectionProfile* profile);

  /// \brief Get the driver to use as specified by this profile.
  ///
  /// It is not required that a profile specify a driver. If the options
  // can be reusable across drivers, then the profile does not need to specify
  /// a driver (if this provides an empty string or nullptr then the driver
  /// must be defined by other means, e.g. by the driver / uri options).
  ///
  /// \param[in] profile The profile to query.
  /// \param[out] driver_name The name of the driver to use, or NULL if not specified.
  /// \param[out] init_func The init function to use for the driver, or NULL if not
  ///   specified.
  /// \param[out] error An optional location to return an error message
  AdbcStatusCode (*GetDriverName)(struct AdbcConnectionProfile* profile,
                                  const char** driver_name, AdbcDriverInitFunc* init_func,
                                  struct AdbcError* error);

  /// \brief Get the string options specified by the profile
  ///
  /// The keys and values returned by this function are owned by the profile
  /// object itself and do not need to be freed or managed by the caller.
  /// They must not be accessed after calling release on the profile.
  ///
  /// The profile can also indicate that a value should be pulled from the environment
  /// by having a value in the form `env_var(ENV_VAR_NAME)`. If the driver
  /// manager encounters a value of this form, it will replace it with the actual value
  /// of the environment variable `ENV_VAR_NAME` before setting the option. This
  /// is only valid for option *values* not *keys*.
  ///
  /// \param[in] profile The profile to query.
  /// \param[out] keys The keys of the options specified by the profile.
  /// \param[out] values The values of the options specified by the profile.
  /// \param[out] num_options The number of options specified by the profile,
  ///   consumers must not access keys or values beyond this count.
  /// \param[out] error An optional location to return an error message
  AdbcStatusCode (*GetOptions)(struct AdbcConnectionProfile* profile, const char*** keys,
                               const char*** values, size_t* num_options,
                               struct AdbcError* error);

  /// \brief Get the integer options specified by the profile
  ///
  /// The keys and values returned by this function are owned by the profile
  /// object itself and do not need to be freed or managed by the caller. They must not be
  /// accessed after calling release on the profile.
  ///
  /// Values returned by this function will be set using the DatabaseSetOptionInt function
  /// on the database object being initialized. If the driver does not support the
  /// DatabaseSetOptionInt function, then options should only be returned as strings.
  ///
  /// \param[in] profile The profile to query.
  /// \param[out] keys The keys of the options specified by the profile.
  /// \param[out] values The values of the options specified by the profile.
  /// \param[out] num_options The number of options specified by the profile,
  ///   consumers must not access keys or values beyond this count.
  /// \param[out] error An optional location to return an error message
  AdbcStatusCode (*GetIntOptions)(struct AdbcConnectionProfile* profile,
                                  const char*** keys, const int64_t** values,
                                  size_t* num_options, struct AdbcError* error);

  /// \brief Get the double options specified by the profile
  ///
  /// The keys and values returned by this function are owned by the profile
  /// object itself and do not need to be freed or managed by the caller. They must not be
  /// accessed after calling release on the profile.
  ///
  /// Values returned by this function will be set using the DatabaseSetOptionDouble
  /// function on the database object being initialized. If the driver does not support
  /// the DatabaseSetOptionDouble function, then options should only be returned as
  /// strings.
  ///
  /// \param[in] profile The profile to query.
  /// \param[out] keys The keys of the options specified by the profile.
  /// \param[out] values The values of the options specified by the profile.
  /// \param[out] num_options The number of options specified by the profile,
  ///   consumers must not access keys or values beyond this count.
  /// \param[out] error An optional location to return an error message
  AdbcStatusCode (*GetDoubleOptions)(struct AdbcConnectionProfile* profile,
                                     const char*** keys, const double** values,
                                     size_t* num_options, struct AdbcError* error);
};

/// \brief Common definition for a connection profile provider
///
/// \param[in] profile_name The name of the profile to load. This is the value of the
///   "profile" option or the profile specified in the URI.
/// \param[in] additional_search_path_list A list of additional paths to search for
///   profiles, delimited by the OS specific path list separator.
/// \param[out] out The profile to return. The caller will take ownership of the profile
///   and is responsible for calling release on it when finished.
/// \param[out] error An optional location to return an error message if necessary.
typedef AdbcStatusCode (*AdbcConnectionProfileProvider)(
    const char* profile_name, const char* additional_search_path_list,
    struct AdbcConnectionProfile* out, struct AdbcError* error);

/// \brief Set a custom connection profile provider for the driver manager.
///
/// If no provider is set, the driver manager will use a default, filesystem-based
/// provider which will look for profiles in the following locations if not given an
/// absolute path to a file:
///
/// 1. The environment variable ADBC_PROFILE_PATH, which is a list of paths to search for
/// profiles.
/// 2. The user-level configuration directory (e.g. ~/.config/adbc/profiles on Linux).
///
/// The filesystem-based profile looks for a file named <profile_name>.toml if there is
/// no extension provided, attempting to parse the toml file for the profile information.
/// If the file is found and parsed successfully, the options specified in the profile
/// which have not already been set will be set as if by DatabaseSetOption just before
/// initialization as part of DatabaseInit.
///
/// For file-based profiles the expected format is as follows:
/// ```toml
/// profile_version = 1
/// driver = "driver_name"
///
/// [Options]
/// option1 = "value1"
/// option2 = 42
/// option3 = 3.14
/// ```
///
/// Boolean options will be converted to string equivalents of "true" or "false".
///
/// \param[in] database The database to set the profile provider for.
/// \param[in] provider The profile provider to use. If NULL, the default filesystem-based
///   provider will be used if a profile is needed.
/// \param[out] error An optional location to return an error message if necessary
ADBC_EXPORT
AdbcStatusCode AdbcDriverManagerDatabaseSetProfileProvider(
    struct AdbcDatabase* database, AdbcConnectionProfileProvider provider,
    struct AdbcError* error);

/// \brief Default Filesystem-based profile provider for the driver manager.
///
/// We expose this so that consumers would be able to write a provider that falls back on
/// the default filesystem-based provider if their custom provider fails to find a
/// profile. This allows for more flexible provider implementations that can still
/// leverage the default behavior when needed.
ADBC_EXPORT
AdbcStatusCode AdbcProfileProviderFilesystem(const char* profile_name,
                                             const char* additional_search_path_list,
                                             struct AdbcConnectionProfile* out,
                                             struct AdbcError* error);

/// @}

#endif  // ADBC_DRIVER_MANAGER_H

#ifdef __cplusplus
}
#endif
