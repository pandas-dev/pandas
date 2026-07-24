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

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>  // Must come first
#endif                // defined(_WIN32)

#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "arrow-adbc/adbc.h"
#include "arrow-adbc/adbc_driver_manager.h"

// Forward declarations and shared utilities for driver manager implementation

// Enums
enum class SearchPathSource {
  kEnv,
  kUser,
  kRegistry,
  kSystem,
  kAdditional,
  kConda,
  kUnset,
  kDoesNotExist,
  kDisabledAtCompileTime,
  kDisabledAtRunTime,
  kOtherError,
};

enum class SearchPathType {
  kManifest,
  kProfile,
};

using SearchPaths = std::vector<std::pair<SearchPathSource, std::filesystem::path>>;

// Structs - forward declarations
struct ParseDriverUriResult {
  std::string_view driver;
  std::optional<std::string_view> uri;
  std::optional<std::string_view> profile;
};

struct DriverInfo {
  std::string manifest_file;
  int64_t manifest_version = 0;
  std::string driver_name;
  std::filesystem::path lib_path;
  std::string entrypoint;
  std::string version;
  std::string source;
};

struct OwnedError {
  struct AdbcError error = ADBC_ERROR_INIT;
  ~OwnedError();
};

inline OwnedError::~OwnedError() {
  if (error.release) {
    error.release(&error);
  }
}

// Platform helpers
#ifdef _WIN32
using char_type = wchar_t;
using string_type = std::wstring;

std::string Utf8Encode(const std::wstring& wstr);
std::wstring Utf8Decode(const std::string& str);
void GetWinError(std::string* buffer);
#else
using char_type = char;
using string_type = std::string;
#endif

struct ManagedLibrary {
  ManagedLibrary() : handle(nullptr) {}
  ManagedLibrary(ManagedLibrary&& other) : handle(other.handle) {
    other.handle = nullptr;
  }
  ManagedLibrary(const ManagedLibrary&) = delete;
  ManagedLibrary& operator=(const ManagedLibrary&) = delete;
  ManagedLibrary& operator=(ManagedLibrary&& other) noexcept {
    this->handle = other.handle;
    other.handle = nullptr;
    return *this;
  }

  ~ManagedLibrary() { Release(); }
  void Release();
  /// \brief Resolve the driver name to a concrete location.
  AdbcStatusCode GetDriverInfo(
      const std::string_view driver_name, const AdbcLoadFlags load_options,
      const std::vector<std::filesystem::path>& additional_search_paths, DriverInfo& info,
      struct AdbcError* error);

  /// \return ADBC_STATUS_NOT_FOUND if the driver shared library could not be
  ///   found (via dlopen) or if a manifest was found but did not contain a
  ///   path for the current platform, ADBC_STATUS_INVALID_ARGUMENT if a
  ///   manifest was found but could not be parsed, ADBC_STATUS_OK otherwise
  ///
  /// May modify search_paths to add error info
  AdbcStatusCode SearchPathsForDriver(const std::filesystem::path& driver_path,
                                      SearchPaths& search_paths, DriverInfo& info,
                                      struct AdbcError* error);

  AdbcStatusCode FindDriver(
      const std::filesystem::path& driver_path, const AdbcLoadFlags load_options,
      const std::vector<std::filesystem::path>& additional_search_paths, DriverInfo& info,
      struct AdbcError* error);

  /// \return ADBC_STATUS_NOT_FOUND if the driver shared library could not be
  ///   found, ADBC_STATUS_OK otherwise
  AdbcStatusCode Load(const string_type& library, const SearchPaths& attempted_paths,
                      struct AdbcError* error);

  AdbcStatusCode Lookup(const char* name, void** func, struct AdbcError* error);

#if defined(_WIN32)
  // The loaded DLL
  HMODULE handle;
#else
  void* handle;
#endif  // defined(_WIN32)
};

// Error handling
void ReleaseError(struct AdbcError* error);
void SetError(struct AdbcError* error, const std::string& message);
void AppendError(struct AdbcError* error, const std::string& message);
void SetError(struct AdbcError* error, struct AdbcError* src_error);

// Utilities
std::string CheckNonPrintableLibraryName(const std::string& name);
bool HasExtension(const std::filesystem::path& path, const std::string& ext);
void AddSearchPathsToError(const SearchPaths& search_paths, const SearchPathType& type,
                           std::string& error_message);
SearchPaths GetEnvPaths(const char_type* env_var);

// Path management
ADBC_EXPORT
std::vector<std::filesystem::path> InternalAdbcParsePath(const std::string_view path);
ADBC_EXPORT
std::filesystem::path InternalAdbcUserConfigDir();
#if !defined(_WIN32)
ADBC_EXPORT
std::filesystem::path InternalAdbcSystemConfigDir();
#endif
ADBC_EXPORT
std::optional<ParseDriverUriResult> InternalAdbcParseDriverUri(std::string_view str);

// Search paths
SearchPaths GetSearchPaths(const AdbcLoadFlags levels);

// Driver loading
AdbcStatusCode LoadDriverManifest(const std::filesystem::path& driver_manifest,
                                  DriverInfo& info, struct AdbcError* error);
ADBC_EXPORT
std::string InternalAdbcDriverManagerDefaultEntrypoint(const std::string& driver);

#ifdef _WIN32
AdbcStatusCode LoadDriverFromRegistry(HKEY root, const std::wstring& driver_name,
                                      DriverInfo& info, struct AdbcError* error);
#endif

// Profile loading
AdbcStatusCode ProcessProfileValue(std::string_view key, std::string_view value,
                                   std::string& out, struct AdbcError* error);

// Initialization
/// Temporary state while the database is being configured.
struct TempDatabase {
  std::unordered_map<std::string, std::string> options;
  std::unordered_map<std::string, std::string> bytes_options;
  std::unordered_map<std::string, int64_t> int_options;
  std::unordered_map<std::string, double> double_options;
  std::string driver;
  std::string entrypoint;
  AdbcDriverInitFunc init_func = nullptr;
  AdbcLoadFlags load_flags = ADBC_LOAD_FLAG_ALLOW_RELATIVE_PATHS;
  std::string additional_manifest_search_path_list;
  std::string additional_profile_search_path_list;
  AdbcConnectionProfileProvider profile_provider = nullptr;
};

// Parse and validate options like uri/profile, with the result being that
// `driver` or `init_func` are populated and options like `profile` are
// removed
ADBC_EXPORT
AdbcStatusCode InternalAdbcParseOptions(TempDatabase* db, struct AdbcError* error);

/// Temporary state while the connection is being configured.
struct TempConnection {
  std::unordered_map<std::string, std::string> options;
  std::unordered_map<std::string, std::string> bytes_options;
  std::unordered_map<std::string, int64_t> int_options;
  std::unordered_map<std::string, double> double_options;
  AdbcConnectionProfile* connection_profile = nullptr;
};

AdbcStatusCode InternalInitializeProfile(TempDatabase* args,
                                         const std::string_view profile,
                                         struct AdbcError* error);

// Driver management
AdbcStatusCode AdbcFindLoadDriver(const char* driver_name, const char* entrypoint,
                                  const int version, const AdbcLoadFlags load_options,
                                  const char* additional_search_path_list,
                                  void* raw_driver, struct AdbcError* error);
// ReleaseDriver is implemented in adbc_driver_manager.cc but used by API implementations
AdbcStatusCode ReleaseDriver(struct AdbcDriver* driver, struct AdbcError* error);

#define INIT_ERROR(ERROR, SOURCE)                                    \
  if ((ERROR) != nullptr &&                                          \
      (ERROR)->vendor_code == ADBC_ERROR_VENDOR_CODE_PRIVATE_DATA) { \
    (ERROR)->private_driver = (SOURCE)->private_driver;              \
  }

#define WRAP_STREAM(EXPR, OUT, SOURCE)                   \
  if (!(OUT)) {                                          \
    /* Happens for ExecuteQuery where out is optional */ \
    return EXPR;                                         \
  }                                                      \
  AdbcStatusCode status_code = EXPR;                     \
  ErrorArrayStreamInit(OUT, (SOURCE)->private_driver);   \
  return status_code;

#define CHECK_STATUS(EXPR)                                \
  if (auto _status = (EXPR); _status != ADBC_STATUS_OK) { \
    return _status;                                       \
  }

#ifdef _WIN32
inline const wchar_t* kAdbcDriverPath = L"ADBC_DRIVER_PATH";
#else
inline const char* kAdbcDriverPath = "ADBC_DRIVER_PATH";
#endif  // _WIN32
