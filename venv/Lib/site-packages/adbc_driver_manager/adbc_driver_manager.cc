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

#if defined(_WIN32)
#include <windows.h>  // Must come first

#ifndef NTDDI_VERSION
#define NTDDI_VERSION 0x0A00000C  // For SHGetKnownFolderPath in ShlObj_core.h in ShlObj.h
#endif

#include <KnownFolders.h>
#include <ShlObj.h>
#include <libloaderapi.h>
#include <string.h>  // _wcsnicmp
#include <strsafe.h>
#include <locale>
#else
#include <dlfcn.h>
#endif  // defined(_WIN32)

#include <toml++/toml.hpp>
#include "adbc_driver_manager_internal.h"
#include "arrow-adbc/adbc.h"
#include "arrow-adbc/adbc_driver_manager.h"
#include "current_arch.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <functional>
#include <regex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std::string_literals;  // NOLINT [build/namespaces]

// NOTE: Error handling, structures, and helper functions are now in internal header
// and shared across all source files. This anonymous namespace only contains
// implementation-specific helpers for this file.

// Generate a note for the error message if the library name has potentially
// non-printable (or really non-ASCII-printable-range) characters.  Oblivious
// to Unicode and locales.
std::string CheckNonPrintableLibraryName(const std::string& name) {
  // We could use std::isprint, but that requires locales; prefer a
  // simpler check for out-of-ASCII-range.
  bool has_non_printable = std::any_of(name.begin(), name.end(), [](char c) {
    int v = static_cast<int>(c);
    return v < 32 || v > 127;
  });
  if (!has_non_printable) return "";

  std::string error_message = "Note: driver name may have non-printable characters: `";
  // TODO(lidavidm): we can simplify with C++20 <format>
  for (char c : name) {
    int v = static_cast<int>(c);
    if (v < 32 || v > 127) {
      error_message += "\\x";
      char buf[3];
      std::snprintf(buf, sizeof(buf), "%02x", v & 0xFF);
      error_message += buf;
    } else {
      error_message += c;
    }
  }
  error_message += "`";
  return error_message;
}

// Platform-specific helpers

#if defined(_WIN32)
/// Append a description of the Windows error to the buffer.
void GetWinError(std::string* buffer) {
  DWORD rc = GetLastError();
  LPVOID message;

  FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                    FORMAT_MESSAGE_IGNORE_INSERTS,
                /*lpSource=*/nullptr, rc, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                reinterpret_cast<LPSTR>(&message), /*nSize=*/0, /*Arguments=*/nullptr);

  (*buffer) += '(';
  (*buffer) += std::to_string(rc);
  (*buffer) += ") ";
  (*buffer) += reinterpret_cast<char*>(message);
  LocalFree(message);
}

#endif  // defined(_WIN32)

// Error handling

#ifdef _WIN32
std::string Utf8Encode(const std::wstring& wstr) {
  if (wstr.empty()) return std::string();
  int size_needed = WideCharToMultiByte(
      CP_UTF8, 0, wstr.data(), static_cast<int>(wstr.size()), NULL, 0, NULL, NULL);
  std::string str_to(size_needed, 0);
  WideCharToMultiByte(CP_UTF8, 0, wstr.data(), static_cast<int>(wstr.size()),
                      str_to.data(), size_needed, NULL, NULL);
  return str_to;
}

std::wstring Utf8Decode(const std::string& str) {
  if (str.empty()) return std::wstring();
  int size_needed =
      MultiByteToWideChar(CP_UTF8, 0, str.data(), static_cast<int>(str.size()), NULL, 0);
  std::wstring wstr_to(size_needed, 0);
  MultiByteToWideChar(CP_UTF8, 0, str.data(), static_cast<int>(str.size()),
                      wstr_to.data(), size_needed);
  return wstr_to;
}
#endif  // _WIN32

/// Hold the driver DLL and the driver release callback in the driver struct.
struct ManagerDriverState {
  // The original release callback
  AdbcStatusCode (*driver_release)(struct AdbcDriver* driver, struct AdbcError* error);

  ManagedLibrary handle;
};

/// Unload the driver DLL.
static AdbcStatusCode ReleaseDriverInternal(struct AdbcDriver* driver,
                                            struct AdbcError* error) {
  AdbcStatusCode status = ADBC_STATUS_OK;

  if (!driver->private_manager) return status;
  ManagerDriverState* state =
      reinterpret_cast<ManagerDriverState*>(driver->private_manager);

  if (state->driver_release) {
    status = state->driver_release(driver, error);
  }
  state->handle.Release();

  driver->private_manager = nullptr;
  delete state;
  return status;
}

// Default stubs

static const char kDefaultEntrypoint[] = "AdbcDriverInit";

// Wrapper to expose ReleaseDriver from anonymous namespace
AdbcStatusCode ReleaseDriver(struct AdbcDriver* driver, struct AdbcError* error) {
  return ReleaseDriverInternal(driver, error);
}

// Utility functions (shared across all source files)
void AddSearchPathsToError(const SearchPaths& search_paths, const SearchPathType& type,
                           std::string& error_message) {
  if (!search_paths.empty()) {
    error_message += "\nAlso searched these paths for";
    if (type == SearchPathType::kManifest) {
      error_message += " manifests:";
    } else if (type == SearchPathType::kProfile) {
      error_message += " profiles:";
    }

    for (const auto& [source, path] : search_paths) {
      error_message += "\n\t";
      switch (source) {
        case SearchPathSource::kEnv:
          if (type == SearchPathType::kManifest) {
            error_message += "ADBC_DRIVER_PATH: ";
          } else if (type == SearchPathType::kProfile) {
            error_message += "ADBC_PROFILE_PATH: ";
          }
          break;
        case SearchPathSource::kUser:
          error_message += "user config dir: ";
          break;
        case SearchPathSource::kRegistry:
          error_message += "Registry: ";
          break;
        case SearchPathSource::kSystem:
          error_message += "system config dir: ";
          break;
        case SearchPathSource::kAdditional:
          error_message += "additional search path: ";
          break;
        case SearchPathSource::kConda:
          error_message += "Conda prefix: ";
          break;
        case SearchPathSource::kUnset:
          error_message += "not set: ";
          break;
        case SearchPathSource::kDoesNotExist:
          error_message += "does not exist: ";
          break;
        case SearchPathSource::kDisabledAtCompileTime:
          error_message += "not enabled at build time: ";
          break;
        case SearchPathSource::kDisabledAtRunTime:
          error_message += "not enabled at run time: ";
          break;
        case SearchPathSource::kOtherError:
          // Don't add any prefix
          break;
      }
      error_message += path.string();
    }
  }
}

// Error handling functions (shared across all source files)
void ReleaseError(struct AdbcError* error) {
  if (error) {
    if (error->message) delete[] error->message;
    error->message = nullptr;
    error->release = nullptr;
  }
}

void SetError(struct AdbcError* error, const std::string& message) {
  static const std::string kPrefix = "[Driver Manager] ";

  if (!error) return;
  if (error->release) error->release(error);

  // Prepend a string to identify driver manager errors
  error->message = new char[kPrefix.size() + message.size() + 1];
  kPrefix.copy(error->message, kPrefix.size());
  message.copy(error->message + kPrefix.size(), message.size());
  error->message[kPrefix.size() + message.size()] = '\0';
  error->release = ReleaseError;
}

void AppendError(struct AdbcError* error, const std::string& message) {
  if (!error) return;
  if (!error->release || !error->message) {
    SetError(error, message);
    return;
  }

  size_t original_length = std::strlen(error->message);
  size_t combined_length = original_length + 1 + message.size() + 1;
  char* new_message = new char[combined_length];
  std::ignore = std::snprintf(new_message, combined_length, "%s\n%s", error->message,
                              message.c_str());

  error->release(error);
  error->message = new_message;
  error->release = ReleaseError;
}

// Copies src_error into error and releases src_error
void SetError(struct AdbcError* error, struct AdbcError* src_error) {
  if (!error) return;
  if (error->release) error->release(error);

  if (src_error->message) {
    size_t message_size = strlen(src_error->message);
    error->message = new char[message_size + 1];  // +1 to include null
    std::memcpy(error->message, src_error->message, message_size);
    error->message[message_size] = '\0';
  } else {
    error->message = nullptr;
  }

  error->release = ReleaseError;
  if (src_error->release) {
    src_error->release(src_error);
  }
}

// Other helpers (intentionally not in an anonymous namespace so they can be tested)
ADBC_EXPORT std::filesystem::path InternalAdbcUserConfigDir() {
  std::filesystem::path config_dir;
#if defined(_WIN32)
  // SHGetFolderPath is just an alias to SHGetKnownFolderPath since Vista
  // so let's just call the updated function.
  PWSTR path = nullptr;
  auto hres = SHGetKnownFolderPath(FOLDERID_LocalAppData, 0, nullptr, &path);
  if (!SUCCEEDED(hres)) {
    return config_dir;
  }

  std::wstring wpath(path);
  std::filesystem::path dir(std::move(wpath));
  if (!dir.empty()) {
    config_dir = std::filesystem::path(dir);
    config_dir /= "ADBC";
  }
#elif defined(__APPLE__)
  auto dir = std::getenv("HOME");
  if (dir) {
    config_dir = std::filesystem::path(dir);
    config_dir /= "Library/Application Support/ADBC";
  }
#elif defined(__linux__)
  auto dir = std::getenv("XDG_CONFIG_HOME");
  if (!dir) {
    dir = std::getenv("HOME");
    if (dir) {
      config_dir = std::filesystem::path(dir) /= ".config";
    }
  } else {
    config_dir = std::filesystem::path(dir);
  }

  if (!config_dir.empty()) {
    config_dir = config_dir / "adbc";
  }
#endif  // defined(_WIN32)

  return config_dir;
}

#if !defined(_WIN32)
std::filesystem::path InternalAdbcSystemConfigDir() {
#if defined(__APPLE__)
  return std::filesystem::path("/Library/Application Support/ADBC/Drivers");
#else
  return std::filesystem::path("/etc/adbc/drivers");
#endif  // defined(__APPLE__)
}
#endif  // !defined(_WIN32)

std::vector<std::filesystem::path> InternalAdbcParsePath(const std::string_view path) {
  std::vector<std::filesystem::path> result;
  if (path.empty()) {
    return result;
  }

#ifdef _WIN32
  constexpr char delimiter = ';';

  // pulling the logic from Go's filepath.SplitList function
  // where windows checks for quoted/escaped sections while splitting
  // but unix doesn't.
  // see
  // https://cs.opensource.google/go/go/+/refs/tags/go1.24.3:src/path/filepath/path_windows.go
  bool in_quotes = false;
  size_t start = 0;
  for (size_t i = 0; i < path.size(); ++i) {
    if (path[i] == '"') {
      in_quotes = !in_quotes;
    } else if (path[i] == delimiter && !in_quotes) {
      result.emplace_back(path.substr(start, i - start));
      start = i + 1;
    }
  }
  result.emplace_back(path.substr(start));
#else
  constexpr char delimiter = ':';

  size_t start = 0;
  size_t end = 0;
  while ((end = path.find(delimiter, start)) != std::string::npos) {
    result.emplace_back(path.substr(start, end - start));
    start = end + 1;
  }
  result.emplace_back(path.substr(start));
#endif  // _WIN32

  // remove empty paths
  result.erase(std::remove_if(result.begin(), result.end(),
                              [](const auto& p) { return p.empty(); }),
               result.end());
  return result;
}

ADBC_EXPORT
std::optional<ParseDriverUriResult> InternalAdbcParseDriverUri(std::string_view str) {
  std::string::size_type pos = str.find(":");
  if (pos == std::string::npos) {
    return std::nullopt;
  }

  std::string_view d = str.substr(0, pos);
  if (str.size() <= pos + 1) {
    return ParseDriverUriResult{d, std::nullopt, std::nullopt};
  }

#ifdef _WIN32
  if (std::filesystem::exists(std::filesystem::path(str))) {
    // No scheme, just a path
    return ParseDriverUriResult{str, std::nullopt, std::nullopt};
  }
#endif

  if (str[pos + 1] == '/') {  // scheme is also driver
    if (d == "profile" && str.size() > pos + 2) {
      // found a profile URI "profile://"
      return ParseDriverUriResult{"", std::nullopt, str.substr(pos + 3)};
    }
    return ParseDriverUriResult{d, str, std::nullopt};
  }

  // driver:scheme:.....
  return ParseDriverUriResult{d, str.substr(pos + 1), std::nullopt};
}

struct ProfileGuard {
  AdbcConnectionProfile profile;
  explicit ProfileGuard() : profile{} {}
  ~ProfileGuard() {
    if (profile.release) {
      profile.release(&profile);
    }
  }
};

ADBC_EXPORT
AdbcStatusCode InternalAdbcParseOptions(TempDatabase* db, struct AdbcError* error) {
  // https://github.com/apache/arrow-adbc/issues/4085

  // init_func always takes precedence
  if (db->init_func) return ADBC_STATUS_OK;

  std::optional<ParseDriverUriResult> parsed_driver;
  std::optional<ParseDriverUriResult> parsed_uri;

  if (!db->driver.empty()) {
    parsed_driver = InternalAdbcParseDriverUri(db->driver);
  }
  if (auto it = db->options.find("uri"); it != db->options.end()) {
    parsed_uri = InternalAdbcParseDriverUri(it->second);
  }

  // do we have `profile`, `profile://` URI, or `profile://` URI in driver?
  {
    const bool have_profile = db->options.find("profile") != db->options.end();
    const bool have_profile_driver = parsed_driver && parsed_driver->profile.has_value();
    const bool have_profile_uri = parsed_uri && parsed_uri->profile.has_value();

    if (static_cast<int>(have_profile) + static_cast<int>(have_profile_uri) +
            static_cast<int>(have_profile_driver) >
        1) {
      SetError(error,
               "Multiple profiles specified; only one of `profile` option, `profile://` "
               "URI, or `profile://` driver is allowed");
      return ADBC_STATUS_INVALID_ARGUMENT;
    }

    // apply profile (may modify options)
    std::string profile;
    if (have_profile) {
      profile = db->options["profile"];
      db->options.erase("profile");
    } else if (have_profile_uri) {
      profile = *parsed_uri->profile;
      db->options.erase("uri");
    } else if (have_profile_driver) {
      profile = *parsed_driver->profile;
      db->driver.clear();
    }

    if (!profile.empty()) {
      auto status = InternalInitializeProfile(db, profile, error);
      if (status != ADBC_STATUS_OK) {
        return status;
      }
    }
  }

  // parse driver for URI. we do this even if it was set via profile
  {
    // Reparse because profile may have modified driver
    std::string owned_driver = db->driver;
    if (auto maybe_driver = InternalAdbcParseDriverUri(owned_driver);
        maybe_driver.has_value()) {
      auto parsed = *maybe_driver;
      // Don't allow recursive profiles (that is the only way we can reach here)
      if (parsed.profile.has_value()) {
        SetError(error, "Profile cannot specify a profile:// URI");
        return ADBC_STATUS_INVALID_ARGUMENT;
      }
      if (parsed.uri.has_value()) {
        db->driver = std::string{parsed.driver};
        // maybe only do this if URI is unset?
        db->options["uri"] = std::string{*parsed.uri};
      }
    }
  }

  // no driver? parse URI
  if (auto it = db->options.find("uri"); db->driver.empty() && it != db->options.end()) {
    // Reparse because profile may have modified uri
    std::string owned_uri = it->second;
    if (auto maybe_uri = InternalAdbcParseDriverUri(it->second); maybe_uri.has_value()) {
      auto parsed = *maybe_uri;
      // Don't allow recursive profiles (that is the only way we can reach here)
      if (parsed.profile.has_value()) {
        SetError(error, "Profile cannot specify a profile:// URI");
        return ADBC_STATUS_INVALID_ARGUMENT;
      }
      if (parsed.uri.has_value()) {
        db->driver = std::string{parsed.driver};
        db->options["uri"] = std::string{*parsed.uri};
      }
    } else if (db->driver.empty()) {
      db->driver = std::move(owned_uri);
    }
  }

  // still no driver? bail
  if (db->driver.empty()) {
    SetError(error, "Must set 'driver' option before AdbcDatabaseInit");
    return ADBC_STATUS_INVALID_ARGUMENT;
  }
  return ADBC_STATUS_OK;
}

AdbcStatusCode InternalInitializeProfile(TempDatabase* args,
                                         const std::string_view profile,
                                         struct AdbcError* error) {
  if (!args->profile_provider) {
    args->profile_provider = AdbcProfileProviderFilesystem;
  }

  ProfileGuard guard{};
  CHECK_STATUS(args->profile_provider(profile.data(),
                                      args->additional_profile_search_path_list.c_str(),
                                      &guard.profile, error));

  const char* driver_name = nullptr;
  AdbcDriverInitFunc init_func = nullptr;
  CHECK_STATUS(
      guard.profile.GetDriverName(&guard.profile, &driver_name, &init_func, error));
  if (driver_name != nullptr && std::strlen(driver_name) > 0) {
    // ensure the parsed driver matches the profile driver if both are specified
    if (!args->driver.empty() && args->driver != driver_name) {
      std::string message = "Profile `";
      message += profile;
      message += "` specifies driver `";
      message += driver_name;
      message += "` which does not match requested driver `" + args->driver + "`";
      SetError(error, message);
      return ADBC_STATUS_INVALID_ARGUMENT;
    }
    args->driver = driver_name;
  }

  if (init_func != nullptr) {
    args->init_func = init_func;
  }

  const char** keys = nullptr;
  const char** values = nullptr;
  size_t num_options = 0;
  const int64_t* int_values = nullptr;
  const double* double_values = nullptr;

  CHECK_STATUS(
      guard.profile.GetOptions(&guard.profile, &keys, &values, &num_options, error));
  for (size_t i = 0; i < num_options; ++i) {
    // use try_emplace so we only add the option if there isn't
    // already an option with the same name
    std::string processed;
    CHECK_STATUS(ProcessProfileValue(keys[i], values[i], processed, error));
    args->options.try_emplace(keys[i], processed);
  }

  CHECK_STATUS(guard.profile.GetIntOptions(&guard.profile, &keys, &int_values,
                                           &num_options, error));
  for (size_t i = 0; i < num_options; ++i) {
    // use try_emplace so we only add the option if there isn't
    // already an option with the same name
    args->int_options.try_emplace(keys[i], int_values[i]);
  }

  CHECK_STATUS(guard.profile.GetDoubleOptions(&guard.profile, &keys, &double_values,
                                              &num_options, error));
  for (size_t i = 0; i < num_options; ++i) {
    // use try_emplace so we only add the option if there isn't already an option with the
    // same name
    args->double_options.try_emplace(keys[i], double_values[i]);
  }

  return ADBC_STATUS_OK;
}

// AdbcStatusCodeMessage moved to adbc_driver_manager_api.cc (API implementation)

AdbcStatusCode AdbcFindLoadDriver(const char* driver_name, const char* entrypoint,
                                  const int version, const AdbcLoadFlags load_options,
                                  const char* additional_search_path_list,
                                  void* raw_driver, struct AdbcError* error) {
  AdbcDriverInitFunc init_func = nullptr;
  std::string error_message;

  switch (version) {
    case ADBC_VERSION_1_0_0:
    case ADBC_VERSION_1_1_0:
      break;
    default:
      SetError(error, "Only ADBC 1.0.0 and 1.1.0 are supported");
      return ADBC_STATUS_NOT_IMPLEMENTED;
  }

  if (!raw_driver) {
    SetError(error, "Driver pointer is null");
    return ADBC_STATUS_INVALID_ARGUMENT;
  }
  if (!driver_name) {
    SetError(error, "Driver name is null");
    return ADBC_STATUS_INVALID_ARGUMENT;
  }

  ManagedLibrary library;
  DriverInfo info;
  if (entrypoint) {
    info.entrypoint = entrypoint;
  }

  std::vector<std::filesystem::path> additional_paths;
  if (additional_search_path_list) {
    additional_paths = InternalAdbcParsePath(additional_search_path_list);
  }

  auto* driver = reinterpret_cast<struct AdbcDriver*>(raw_driver);

  AdbcStatusCode status =
      library.GetDriverInfo(driver_name, load_options, additional_paths, info, error);
  if (status != ADBC_STATUS_OK) {
    driver->release = nullptr;
    return status;
  }

  void* load_handle = nullptr;
  if (!info.entrypoint.empty()) {
    status = library.Lookup(info.entrypoint.c_str(), &load_handle, error);
  } else {
    auto name = InternalAdbcDriverManagerDefaultEntrypoint(info.lib_path.string());
    assert(!name.empty());
    status = library.Lookup(name.c_str(), &load_handle, error);
    if (status != ADBC_STATUS_OK) {
      status = library.Lookup(kDefaultEntrypoint, &load_handle, error);
    }
  }

  if (status != ADBC_STATUS_OK) {
    library.Release();
    return status;
  }
  init_func = reinterpret_cast<AdbcDriverInitFunc>(load_handle);

  status = AdbcLoadDriverFromInitFunc(init_func, version, driver, error);
  if (status == ADBC_STATUS_OK) {
    ManagerDriverState* state = new ManagerDriverState;
    state->driver_release = driver->release;
    state->handle = std::move(library);
    driver->release = &ReleaseDriver;
    driver->private_manager = state;
  } else {
    library.Release();
  }
  return status;
}
