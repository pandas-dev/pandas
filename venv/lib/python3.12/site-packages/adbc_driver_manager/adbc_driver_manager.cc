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
#define NOMINMAX
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
#include "arrow-adbc/adbc.h"
#include "arrow-adbc/adbc_driver_manager.h"
#include "current_arch.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std::string_literals;  // NOLINT [build/namespaces]

ADBC_EXPORT
std::vector<std::filesystem::path> InternalAdbcParsePath(const std::string_view path);
ADBC_EXPORT
std::filesystem::path InternalAdbcUserConfigDir();

namespace {

/// \brief Where a search path came from (for error reporting)
enum class SearchPathSource {
  kEnv,
  kUser,
  kRegistry,
  kSystem,
  kAdditional,
  kConda,
  kUnset,
  kDoesNotExist,
  kDisabled,
  kOtherError,
};

using SearchPaths = std::vector<std::pair<SearchPathSource, std::filesystem::path>>;

void AddSearchPathsToError(const SearchPaths& search_paths, std::string& error_message) {
  if (!search_paths.empty()) {
    error_message += "\nAlso searched these paths for manifests:";
    for (const auto& [source, path] : search_paths) {
      error_message += "\n\t";
      switch (source) {
        case SearchPathSource::kEnv:
          error_message += "ADBC_DRIVER_PATH: ";
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
        case SearchPathSource::kDisabled:
          error_message += "not enabled at build time: ";
          break;
        case SearchPathSource::kOtherError:
          // Don't add any prefix
          break;
      }
      error_message += path.string();
    }
  }
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

// Copies src_error into error and releases src_error
void SetError(struct AdbcError* error, struct AdbcError* src_error) {
  if (!error) return;
  if (error->release) error->release(error);

  if (src_error->message) {
    size_t message_size = strlen(src_error->message);
    error->message = new char[message_size];
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

struct OwnedError {
  struct AdbcError error = ADBC_ERROR_INIT;

  ~OwnedError() {
    if (error.release) {
      error.release(&error);
    }
  }
};

#ifdef _WIN32
using char_type = wchar_t;

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

#else
using char_type = char;
#endif  // _WIN32

/// \brief The location and entrypoint of a resolved driver.
struct DriverInfo {
  std::string manifest_file;
  int64_t manifest_version;
  std::string driver_name;
  std::filesystem::path lib_path;
  std::string entrypoint;

  std::string version;
  std::string source;
};

#ifdef _WIN32
class RegistryKey {
 public:
  RegistryKey(HKEY root, const std::wstring_view subkey) noexcept
      : root_(root), key_(nullptr) {
    status_ = RegOpenKeyExW(root_, subkey.data(), 0, KEY_READ, &key_);
  }

  ~RegistryKey() {
    if (is_open() && key_ != nullptr) {
      RegCloseKey(key_);
      key_ = nullptr;
      status_ = ERROR_REGISTRY_IO_FAILED;
    }
  }

  HKEY key() const { return key_; }
  bool is_open() const { return status_ == ERROR_SUCCESS; }
  LSTATUS status() const { return status_; }

  std::wstring GetString(const std::wstring& name, std::wstring default_value) {
    if (!is_open()) return default_value;

    DWORD type = REG_SZ;
    DWORD size = 0;
    auto result = RegQueryValueExW(key_, name.data(), nullptr, &type, nullptr, &size);
    if (result != ERROR_SUCCESS) return default_value;
    if (type != REG_SZ) return default_value;

    std::wstring value(size, '\0');
    result = RegQueryValueExW(key_, name.data(), nullptr, &type,
                              reinterpret_cast<LPBYTE>(value.data()), &size);
    if (result != ERROR_SUCCESS) return default_value;
    return value;
  }

  int32_t GetInt(const std::wstring& name, const int32_t default_value) {
    if (!is_open()) return default_value;

    DWORD dwValue;
    DWORD dataSize = sizeof(dwValue);
    DWORD valueType;
    auto result = RegQueryValueExW(key_, name.data(), nullptr, &valueType,
                                   (LPBYTE)&dwValue, &dataSize);
    if (result != ERROR_SUCCESS) return default_value;
    if (valueType != REG_DWORD) return default_value;
    return static_cast<int32_t>(dwValue);
  }

 private:
  HKEY root_;
  HKEY key_;
  LSTATUS status_;
};

AdbcStatusCode LoadDriverFromRegistry(HKEY root, const std::wstring& driver_name,
                                      DriverInfo& info, struct AdbcError* error) {
  // N.B. start all error messages with the subkey so that the calling code
  // can prepend the name of 'root' to the error message (easier than trying
  // to invoke win32 API to get the name of the HKEY)
  static const LPCWSTR kAdbcDriverRegistry = L"SOFTWARE\\ADBC\\Drivers";
  RegistryKey drivers_key(root, kAdbcDriverRegistry);
  if (!drivers_key.is_open()) {
    std::string error_message = "SOFTWARE\\ADBC\\DRIVERS not found"s;
    SetError(error, std::move(error_message));
    return ADBC_STATUS_NOT_FOUND;
  }

  RegistryKey dkey(drivers_key.key(), driver_name);
  if (!dkey.is_open()) {
    std::string error_message = "SOFTWARE\\ADBC\\DRIVERS has no entry for driver \""s;
    error_message += Utf8Encode(driver_name);
    error_message += "\""s;
    SetError(error, std::move(error_message));
    return ADBC_STATUS_NOT_FOUND;
  }

  info.driver_name = Utf8Encode(dkey.GetString(L"name", L""));
  info.manifest_version = int64_t(dkey.GetInt(L"manifest_version", 1));
  if (info.manifest_version != 1) {
    SetError(error, "Driver manifest version '" + std::to_string(info.manifest_version) +
                        "' is not supported by this driver manager.");
    return ADBC_STATUS_INVALID_ARGUMENT;
  }

  info.entrypoint = Utf8Encode(dkey.GetString(L"entrypoint", L""));
  info.version = Utf8Encode(dkey.GetString(L"version", L""));
  info.source = Utf8Encode(dkey.GetString(L"source", L""));
  info.lib_path = std::filesystem::path(dkey.GetString(L"driver", L""));
  if (info.lib_path.empty()) {
    std::string error_message = "SOFTWARE\\ADBC\\DRIVERS\\"s;
    error_message += Utf8Encode(driver_name);
    error_message += " has no driver path"s;
    SetError(error, std::move(error_message));
    return ADBC_STATUS_NOT_FOUND;
  }
  return ADBC_STATUS_OK;
}
#endif  // _WIN32

/// \return ADBC_STATUS_NOT_FOUND if the manifest does not contain a driver
///   path for this platform, ADBC_STATUS_INVALID_ARGUMENT if the manifest
///   could not be parsed, ADBC_STATUS_OK otherwise (`info` will be populated)
AdbcStatusCode LoadDriverManifest(const std::filesystem::path& driver_manifest,
                                  DriverInfo& info, struct AdbcError* error) {
  toml::table config;
  try {
    config = toml::parse_file(driver_manifest.native());
  } catch (const toml::parse_error& err) {
    // Despite the name, this exception covers IO errors too.  Hence, we can't
    // differentiate between bad syntax and other I/O error.
    std::string message = "Could not open manifest. ";
    message += err.what();
    message += ". Manifest: ";
    message += driver_manifest.string();
    SetError(error, std::move(message));
    return ADBC_STATUS_INVALID_ARGUMENT;
  }

  info.manifest_file = driver_manifest.string();
  info.driver_name = config["name"].value_or(""s);
  info.manifest_version = config["manifest_version"].value_or(int64_t(1));
  if (info.manifest_version != 1) {
    SetError(error, "Driver manifest version '" + std::to_string(info.manifest_version) +
                        "' is not supported by this driver manager.");
    return ADBC_STATUS_INVALID_ARGUMENT;
  }

  info.entrypoint = config.at_path("Driver.entrypoint").value_or(""s);
  info.version = config["version"].value_or(""s);
  info.source = config["source"].value_or(""s);

  auto entrypoint = config.at_path("Driver.entrypoint");
  if (entrypoint) {
    if (auto* ep = entrypoint.as_string()) {
      info.entrypoint = ep->get();
    } else {
      SetError(error, "Driver entrypoint not a string in manifest '"s +
                          driver_manifest.string() + "'"s);
      return ADBC_STATUS_INVALID_ARGUMENT;
    }
  }

  auto driver = config.at_path("Driver.shared");
  if (toml::table* platforms = driver.as_table()) {
    auto view = platforms->at_path(adbc::CurrentArch());
    if (!view) {
      std::string message = "Driver path not found in manifest '";
      message += driver_manifest.string();
      message += "' for current architecture '";
      message += adbc::CurrentArch();
      message += "'. Architectures found:";
      for (const auto& [key, val] : *platforms) {
        message += " ";
        message += key;
      }
      SetError(error, std::move(message));
      return ADBC_STATUS_NOT_FOUND;
    } else if (auto* path = view.as_string()) {
      if (path->get().empty()) {
        std::string message = "Driver path is an empty string in manifest '";
        message += driver_manifest.string();
        message += "' for current architecture '";
        message += adbc::CurrentArch();
        message += "'";
        SetError(error, std::move(message));
        return ADBC_STATUS_INVALID_ARGUMENT;
      }

      info.lib_path = path->get();
      return ADBC_STATUS_OK;
    } else {
      std::string message = "Driver path not found in manifest '";
      message += driver_manifest.string();
      message += "' for current architecture '";
      message += adbc::CurrentArch();
      message += "'. Value was not a string";
      SetError(error, std::move(message));
      return ADBC_STATUS_INVALID_ARGUMENT;
    }
    return ADBC_STATUS_OK;
  } else if (auto* path = driver.as_string()) {
    info.lib_path = path->get();
    if (info.lib_path.empty()) {
      SetError(error, "Driver path is an empty string in manifest '"s +
                          driver_manifest.string() + "'"s);
      return ADBC_STATUS_INVALID_ARGUMENT;
    }
    return ADBC_STATUS_OK;
  }
  SetError(error, "Driver path not defined in manifest '"s + driver_manifest.string() +
                      "'. `Driver.shared` must be a string or table"s);
  return ADBC_STATUS_INVALID_ARGUMENT;
}

SearchPaths GetEnvPaths(const char_type* env_var) {
#ifdef _WIN32
  size_t required_size;

  _wgetenv_s(&required_size, NULL, 0, env_var);
  if (required_size == 0) {
    return {};
  }

  std::wstring path_var;
  path_var.resize(required_size);
  _wgetenv_s(&required_size, path_var.data(), required_size, env_var);
  // Remove null terminator
  path_var.resize(required_size - 1);
  auto path = Utf8Encode(path_var);
#else
  const char* path_var = std::getenv(env_var);
  if (!path_var) {
    return {};
  }
  std::string path(path_var);
#endif  // _WIN32
  SearchPaths paths;
  for (auto path : InternalAdbcParsePath(path)) {
    paths.emplace_back(SearchPathSource::kEnv, path);
  }
  return paths;
}

#ifdef _WIN32
static const wchar_t* kAdbcDriverPath = L"ADBC_DRIVER_PATH";
#else
static const char* kAdbcDriverPath = "ADBC_DRIVER_PATH";
#endif  // _WIN32

SearchPaths GetSearchPaths(const AdbcLoadFlags levels) {
  SearchPaths paths;
  if (levels & ADBC_LOAD_FLAG_SEARCH_ENV) {
    // Check the ADBC_DRIVER_PATH environment variable
    paths = GetEnvPaths(kAdbcDriverPath);
  }

  if (levels & ADBC_LOAD_FLAG_SEARCH_USER) {
    // Check the user configuration directory
    std::filesystem::path user_config_dir = InternalAdbcUserConfigDir();
    if (!user_config_dir.empty() && std::filesystem::exists(user_config_dir)) {
      paths.emplace_back(SearchPathSource::kUser, std::move(user_config_dir));
    } else {
      paths.emplace_back(SearchPathSource::kDoesNotExist, std::move(user_config_dir));
    }
  }

  if (levels & ADBC_LOAD_FLAG_SEARCH_SYSTEM) {
    // System level behavior for Windows is to search the registry keys so we
    // only need to check for macOS and fall back to Unix-like behavior as long
    // as we're not on Windows
#if defined(__APPLE__)
    const std::filesystem::path system_config_dir(
        "/Library/Application Support/ADBC/Drivers");
    if (std::filesystem::exists(system_config_dir)) {
      paths.emplace_back(SearchPathSource::kSystem, std::move(system_config_dir));
    } else {
      paths.emplace_back(SearchPathSource::kDoesNotExist, std::move(system_config_dir));
    }
#elif !defined(_WIN32)
    const std::filesystem::path system_config_dir("/etc/adbc/drivers");
    if (std::filesystem::exists(system_config_dir)) {
      paths.emplace_back(SearchPathSource::kSystem, std::move(system_config_dir));
    } else {
      paths.emplace_back(SearchPathSource::kDoesNotExist, std::move(system_config_dir));
    }
#endif  // defined(__APPLE__)
  }

  return paths;
}

bool HasExtension(const std::filesystem::path& path, const std::string& ext) {
#ifdef _WIN32
  auto wext = Utf8Decode(ext);
  auto path_ext = path.extension().native();
  return path_ext.size() == wext.size() &&
         _wcsnicmp(path_ext.data(), wext.data(), wext.size()) == 0;
#else
  return path.extension() == ext;
#endif  // _WIN32
}

/// A driver DLL.
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

  void Release() {
    // TODO(apache/arrow-adbc#204): causes tests to segfault.  Need to
    // refcount the driver DLL; also, errors may retain a reference to
    // release() from the DLL - how to handle this?  It's unlikely we can
    // actually do this - in general shared libraries are not safe to unload.
  }

  /// \brief Resolve the driver name to a concrete location.
  AdbcStatusCode GetDriverInfo(
      const std::string_view driver_name, const AdbcLoadFlags load_options,
      const std::vector<std::filesystem::path>& additional_search_paths, DriverInfo& info,
      struct AdbcError* error) {
    if (driver_name.empty()) {
      SetError(error, "Driver name is empty");
      return ADBC_STATUS_INVALID_ARGUMENT;
    }

    // First try to treat the given driver name as a path to a manifest or shared library
    std::filesystem::path driver_path(driver_name);
    const bool allow_relative_paths = load_options & ADBC_LOAD_FLAG_ALLOW_RELATIVE_PATHS;
    if (driver_path.has_extension()) {
      if (driver_path.is_relative() && !allow_relative_paths) {
        SetError(error, "Driver path is relative and relative paths are not allowed");
        return ADBC_STATUS_INVALID_ARGUMENT;
      }

      if (HasExtension(driver_path, ".toml")) {
        // if the extension is .toml, attempt to load the manifest
        // erroring if we fail

        auto status = LoadDriverManifest(driver_path, info, error);
        if (status == ADBC_STATUS_OK) {
          return Load(info.lib_path.c_str(), {}, error);
        }
        return status;
      }

      // if the extension is not .toml, then just try to load the provided
      // path as if it was an absolute path to a driver library
      return Load(driver_path.c_str(), {}, error);
    }

    if (driver_path.is_absolute()) {
      // if we have an absolute path without an extension, first see if there's a
      // toml file with the same name.
      driver_path.replace_extension(".toml");
      if (std::filesystem::exists(driver_path)) {
        auto status = LoadDriverManifest(driver_path, info, error);
        if (status == ADBC_STATUS_OK) {
          return Load(info.lib_path.c_str(), {}, error);
        }
      }

      driver_path.replace_extension("");
      info.lib_path = driver_path;
      // otherwise just try to load the provided path as if it was an absolute path
      return Load(driver_path.c_str(), {}, error);
    }

    if (driver_path.has_extension()) {
      if (driver_path.is_relative() && !allow_relative_paths) {
        SetError(error, "Driver path is relative and relative paths are not allowed");
        return ADBC_STATUS_INVALID_ARGUMENT;
      }

#if defined(_WIN32)
      static const std::string kPlatformLibrarySuffix = ".dll";
#elif defined(__APPLE__)
      static const std::string kPlatformLibrarySuffix = ".dylib";
#else
      static const std::string kPlatformLibrarySuffix = ".so";
#endif  // defined(_WIN32)
      if (HasExtension(driver_path, kPlatformLibrarySuffix)) {
        info.lib_path = driver_path;
        return Load(driver_path.c_str(), {}, error);
      }

      SetError(error, "Driver name has unrecognized extension: " +
                          driver_path.extension().string());
      return ADBC_STATUS_INVALID_ARGUMENT;
    }

    // not an absolute path, no extension. Let's search the configured paths
    // based on the options
    return FindDriver(driver_path, load_options, additional_search_paths, info, error);
  }

  /// \return ADBC_STATUS_NOT_FOUND if the driver shared library could not be
  ///   found (via dlopen) or if a manifest was found but did not contain a
  ///   path for the current platform, ADBC_STATUS_INVALID_ARGUMENT if a
  ///   manifest was found but could not be parsed, ADBC_STATUS_OK otherwise
  ///
  /// May modify search_paths to add error info
  AdbcStatusCode SearchPathsForDriver(const std::filesystem::path& driver_path,
                                      SearchPaths& search_paths, DriverInfo& info,
                                      struct AdbcError* error) {
    SearchPaths extra_debug_info;
    for (const auto& [source, search_path] : search_paths) {
      if (source == SearchPathSource::kRegistry || source == SearchPathSource::kUnset ||
          source == SearchPathSource::kDoesNotExist ||
          source == SearchPathSource::kDisabled ||
          source == SearchPathSource::kOtherError) {
        continue;
      }
      std::filesystem::path full_path = search_path / driver_path;

      // check for toml first, then dll
      full_path.replace_extension(".toml");
      if (std::filesystem::exists(full_path)) {
        OwnedError intermediate_error;

        auto status = LoadDriverManifest(full_path, info, &intermediate_error.error);
        if (status == ADBC_STATUS_OK) {
          // Don't pass attempted_paths here; we'll generate the error at a higher level
          status = Load(info.lib_path.c_str(), {}, &intermediate_error.error);
          if (status == ADBC_STATUS_OK) {
            return status;
          }
          std::string message = "found ";
          message += full_path.string();
          if (intermediate_error.error.message) {
            message += " but: ";
            message += intermediate_error.error.message;
          } else {
            message += " could not load the driver it specified";
          }
          extra_debug_info.emplace_back(SearchPathSource::kOtherError,
                                        std::move(message));
          search_paths.insert(search_paths.end(), extra_debug_info.begin(),
                              extra_debug_info.end());
          return status;
        } else if (status == ADBC_STATUS_INVALID_ARGUMENT) {
          // The manifest was invalid. Don't ignore that!
          search_paths.insert(search_paths.end(), extra_debug_info.begin(),
                              extra_debug_info.end());
          if (intermediate_error.error.message) {
            std::string error_message = intermediate_error.error.message;
            AddSearchPathsToError(search_paths, error_message);
            SetError(error, std::move(error_message));
          }
          return status;
        }
        // Should be NOT_FOUND otherwise
        std::string message = "found ";
        message += full_path.string();
        if (intermediate_error.error.message) {
          message += " but: ";
          message += intermediate_error.error.message;
        } else {
          message += " which did not define a driver for this platform";
        }

        extra_debug_info.emplace_back(SearchPathSource::kOtherError, std::move(message));
      }

      // remove the .toml extension; Load will add the DLL/SO/DYLIB suffix
      full_path.replace_extension("");
      // Don't pass error here - it'll be suppressed anyways
      auto status = Load(full_path.c_str(), {}, nullptr);
      if (status == ADBC_STATUS_OK) {
        return status;
      }
    }

    search_paths.insert(search_paths.end(), extra_debug_info.begin(),
                        extra_debug_info.end());
    return ADBC_STATUS_NOT_FOUND;
  }

  AdbcStatusCode FindDriver(
      const std::filesystem::path& driver_path, const AdbcLoadFlags load_options,
      const std::vector<std::filesystem::path>& additional_search_paths, DriverInfo& info,
      struct AdbcError* error) {
    if (driver_path.empty()) {
      SetError(error, "Driver path is empty");
      return ADBC_STATUS_INVALID_ARGUMENT;
    }

    SearchPaths search_paths;
    {
      // First search the paths in the env var `ADBC_DRIVER_PATH`.
      // Then search the runtime application-defined additional search paths.
      search_paths = GetSearchPaths(load_options & ADBC_LOAD_FLAG_SEARCH_ENV);
      if (search_paths.empty()) {
        search_paths.emplace_back(SearchPathSource::kUnset, "ADBC_DRIVER_PATH");
      }
      for (const auto& path : additional_search_paths) {
        search_paths.emplace_back(SearchPathSource::kAdditional, path);
      }

#if ADBC_CONDA_BUILD
      // Then, if this is a conda build, search in the conda environment if
      // it is activated.
      if (load_options & ADBC_LOAD_FLAG_SEARCH_ENV) {
#ifdef _WIN32
        const wchar_t* conda_name = L"CONDA_PREFIX";
#else
        const char* conda_name = "CONDA_PREFIX";
#endif  // _WIN32
        auto venv = GetEnvPaths(conda_name);
        if (!venv.empty()) {
          for (const auto& venv_path : venv) {
            search_paths.emplace_back(SearchPathSource::kConda,
                                      venv_path / "etc" / "adbc" / "drivers");
          }
        }
      }
#else
      if (load_options & ADBC_LOAD_FLAG_SEARCH_ENV) {
        search_paths.emplace_back(SearchPathSource::kDisabled, "Conda prefix");
      }
#endif  // ADBC_CONDA_BUILD

      auto status = SearchPathsForDriver(driver_path, search_paths, info, error);
      if (status != ADBC_STATUS_NOT_FOUND) {
        // If NOT_FOUND, then keep searching; if OK or INVALID_ARGUMENT, stop
        return status;
      }
    }

    // We searched environment paths and additional search paths (if they
    // exist), so now search the rest.
#ifdef _WIN32
    // On Windows, check registry keys, not just search paths.
    if (load_options & ADBC_LOAD_FLAG_SEARCH_USER) {
      // Check the user registry for the driver.
      auto status =
          LoadDriverFromRegistry(HKEY_CURRENT_USER, driver_path.native(), info, error);
      if (status == ADBC_STATUS_OK) {
        return Load(info.lib_path.c_str(), {}, error);
      }
      if (error && error->message) {
        std::string message = "HKEY_CURRENT_USER\\"s;
        message += error->message;
        search_paths.emplace_back(SearchPathSource::kRegistry, std::move(message));
      } else {
        search_paths.emplace_back(SearchPathSource::kRegistry,
                                  "not found in HKEY_CURRENT_USER");
      }

      auto user_paths = GetSearchPaths(ADBC_LOAD_FLAG_SEARCH_USER);
      status = SearchPathsForDriver(driver_path, user_paths, info, error);
      if (status != ADBC_STATUS_NOT_FOUND) {
        return status;
      }
      search_paths.insert(search_paths.end(), user_paths.begin(), user_paths.end());
    }

    if (load_options & ADBC_LOAD_FLAG_SEARCH_SYSTEM) {
      // Check the system registry for the driver.
      auto status =
          LoadDriverFromRegistry(HKEY_LOCAL_MACHINE, driver_path.native(), info, error);
      if (status == ADBC_STATUS_OK) {
        return Load(info.lib_path.c_str(), {}, error);
      }
      if (error && error->message) {
        std::string message = "HKEY_LOCAL_MACHINE\\"s;
        message += error->message;
        search_paths.emplace_back(SearchPathSource::kRegistry, std::move(message));
      } else {
        search_paths.emplace_back(SearchPathSource::kRegistry,
                                  "not found in HKEY_LOCAL_MACHINE");
      }

      auto system_paths = GetSearchPaths(ADBC_LOAD_FLAG_SEARCH_SYSTEM);
      status = SearchPathsForDriver(driver_path, system_paths, info, error);
      if (status != ADBC_STATUS_NOT_FOUND) {
        return status;
      }
      search_paths.insert(search_paths.end(), system_paths.begin(), system_paths.end());
    }

    info.lib_path = driver_path;
    return Load(driver_path.c_str(), search_paths, error);
#else
    // Otherwise, search the configured paths.
    SearchPaths more_search_paths =
        GetSearchPaths(load_options & ~ADBC_LOAD_FLAG_SEARCH_ENV);
    auto status = SearchPathsForDriver(driver_path, more_search_paths, info, error);
    if (status == ADBC_STATUS_NOT_FOUND) {
      // If we reach here, we didn't find the driver in any of the paths
      // so let's just attempt to load it as default behavior.
      search_paths.insert(search_paths.end(), more_search_paths.begin(),
                          more_search_paths.end());
      info.lib_path = driver_path;
      return Load(driver_path.c_str(), search_paths, error);
    }
    return status;
#endif  // _WIN32
  }

  /// \return ADBC_STATUS_NOT_FOUND if the driver shared library could not be
  ///   found, ADBC_STATUS_OK otherwise
  AdbcStatusCode Load(const char_type* library, const SearchPaths& attempted_paths,
                      struct AdbcError* error) {
    std::string error_message;
#if defined(_WIN32)
    HMODULE handle = LoadLibraryExW(library, NULL, 0);
    if (!handle) {
      error_message += Utf8Encode(library);
      error_message += ": LoadLibraryExW() failed: ";
      GetWinError(&error_message);

      std::wstring full_driver_name = library;
      full_driver_name += L".dll";
      handle = LoadLibraryExW(full_driver_name.c_str(), NULL, 0);
      if (!handle) {
        error_message += '\n';
        error_message += Utf8Encode(full_driver_name);
        error_message += ": LoadLibraryExW() failed: ";
        GetWinError(&error_message);
      }
    }
    if (!handle) {
      AddSearchPathsToError(attempted_paths, error_message);
      SetError(error, error_message);
      return ADBC_STATUS_NOT_FOUND;
    } else {
      this->handle = handle;
    }
#else
    static const std::string kPlatformLibraryPrefix = "lib";
#if defined(__APPLE__)
    static const std::string kPlatformLibrarySuffix = ".dylib";
#else
    static const std::string kPlatformLibrarySuffix = ".so";
#endif  // defined(__APPLE__)

    void* handle = dlopen(library, RTLD_NOW | RTLD_LOCAL);
    if (!handle) {
      error_message = "dlopen() failed: ";
      error_message += dlerror();

      // If applicable, append the shared library prefix/extension and
      // try again (this way you don't have to hardcode driver names by
      // platform in the application)
      const std::string driver_str = library;

      std::string full_driver_name;
      if (driver_str.size() < kPlatformLibraryPrefix.size() ||
          driver_str.compare(0, kPlatformLibraryPrefix.size(), kPlatformLibraryPrefix) !=
              0) {
        full_driver_name += kPlatformLibraryPrefix;
      }
      full_driver_name += library;
      if (driver_str.size() < kPlatformLibrarySuffix.size() ||
          driver_str.compare(full_driver_name.size() - kPlatformLibrarySuffix.size(),
                             kPlatformLibrarySuffix.size(),
                             kPlatformLibrarySuffix) != 0) {
        full_driver_name += kPlatformLibrarySuffix;
      }
      handle = dlopen(full_driver_name.c_str(), RTLD_NOW | RTLD_LOCAL);
      if (!handle) {
        error_message += "\ndlopen() failed: ";
        error_message += dlerror();
      }
    }
    if (handle) {
      this->handle = handle;
    } else {
      AddSearchPathsToError(attempted_paths, error_message);
      SetError(error, error_message);
      return ADBC_STATUS_NOT_FOUND;
    }
#endif  // defined(_WIN32)
    return ADBC_STATUS_OK;
  }

  AdbcStatusCode Lookup(const char* name, void** func, struct AdbcError* error) {
#if defined(_WIN32)
    void* load_handle = reinterpret_cast<void*>(GetProcAddress(handle, name));
    if (!load_handle) {
      std::string message = "GetProcAddress(";
      message += name;
      message += ") failed: ";
      GetWinError(&message);
      SetError(error, message);
      return ADBC_STATUS_INTERNAL;
    }
#else
    void* load_handle = dlsym(handle, name);
    if (!load_handle) {
      std::string message = "dlsym(";
      message += name;
      message += ") failed: ";
      message += dlerror();
      SetError(error, message);
      return ADBC_STATUS_INTERNAL;
    }
#endif  // defined(_WIN32)
    *func = load_handle;
    return ADBC_STATUS_OK;
  }

#if defined(_WIN32)
  // The loaded DLL
  HMODULE handle;
#else
  void* handle;
#endif  // defined(_WIN32)
};

/// Hold the driver DLL and the driver release callback in the driver struct.
struct ManagerDriverState {
  // The original release callback
  AdbcStatusCode (*driver_release)(struct AdbcDriver* driver, struct AdbcError* error);

  ManagedLibrary handle;
};

/// Unload the driver DLL.
static AdbcStatusCode ReleaseDriver(struct AdbcDriver* driver, struct AdbcError* error) {
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

// ArrowArrayStream wrapper to support AdbcErrorFromArrayStream

struct ErrorArrayStream {
  struct ArrowArrayStream stream;
  struct AdbcDriver* private_driver;
};

void ErrorArrayStreamRelease(struct ArrowArrayStream* stream) {
  if (stream->release != ErrorArrayStreamRelease || !stream->private_data) return;

  auto* private_data = reinterpret_cast<struct ErrorArrayStream*>(stream->private_data);
  private_data->stream.release(&private_data->stream);
  delete private_data;
  std::memset(stream, 0, sizeof(*stream));
}

const char* ErrorArrayStreamGetLastError(struct ArrowArrayStream* stream) {
  if (stream->release != ErrorArrayStreamRelease || !stream->private_data) return nullptr;
  auto* private_data = reinterpret_cast<struct ErrorArrayStream*>(stream->private_data);
  return private_data->stream.get_last_error(&private_data->stream);
}

int ErrorArrayStreamGetNext(struct ArrowArrayStream* stream, struct ArrowArray* array) {
  if (stream->release != ErrorArrayStreamRelease || !stream->private_data) return EINVAL;
  auto* private_data = reinterpret_cast<struct ErrorArrayStream*>(stream->private_data);
  return private_data->stream.get_next(&private_data->stream, array);
}

int ErrorArrayStreamGetSchema(struct ArrowArrayStream* stream,
                              struct ArrowSchema* schema) {
  if (stream->release != ErrorArrayStreamRelease || !stream->private_data) return EINVAL;
  auto* private_data = reinterpret_cast<struct ErrorArrayStream*>(stream->private_data);
  return private_data->stream.get_schema(&private_data->stream, schema);
}

// Default stubs

int ErrorGetDetailCount(const struct AdbcError* error) { return 0; }

struct AdbcErrorDetail ErrorGetDetail(const struct AdbcError* error, int index) {
  return {nullptr, nullptr, 0};
}

const struct AdbcError* ErrorFromArrayStream(struct ArrowArrayStream* stream,
                                             AdbcStatusCode* status) {
  return nullptr;
}

void ErrorArrayStreamInit(struct ArrowArrayStream* out,
                          struct AdbcDriver* private_driver) {
  if (!out || !out->release ||
      // Don't bother wrapping if driver didn't claim support
      private_driver->ErrorFromArrayStream == ErrorFromArrayStream) {
    return;
  }
  struct ErrorArrayStream* private_data = new ErrorArrayStream;
  private_data->stream = *out;
  private_data->private_driver = private_driver;
  out->get_last_error = ErrorArrayStreamGetLastError;
  out->get_next = ErrorArrayStreamGetNext;
  out->get_schema = ErrorArrayStreamGetSchema;
  out->release = ErrorArrayStreamRelease;
  out->private_data = private_data;
}

AdbcStatusCode DatabaseGetOption(struct AdbcDatabase* database, const char* key,
                                 char* value, size_t* length, struct AdbcError* error) {
  SetError(error, "AdbcDatabaseGetOption not implemented");
  return ADBC_STATUS_NOT_FOUND;
}

AdbcStatusCode DatabaseGetOptionBytes(struct AdbcDatabase* database, const char* key,
                                      uint8_t* value, size_t* length,
                                      struct AdbcError* error) {
  SetError(error, "AdbcDatabaseGetOptionBytes not implemented");
  return ADBC_STATUS_NOT_FOUND;
}

AdbcStatusCode DatabaseGetOptionInt(struct AdbcDatabase* database, const char* key,
                                    int64_t* value, struct AdbcError* error) {
  SetError(error, "AdbcDatabaseGetOptionInt not implemented");
  return ADBC_STATUS_NOT_FOUND;
}

AdbcStatusCode DatabaseGetOptionDouble(struct AdbcDatabase* database, const char* key,
                                       double* value, struct AdbcError* error) {
  SetError(error, "AdbcDatabaseGetOptionDouble not implemented");
  return ADBC_STATUS_NOT_FOUND;
}

AdbcStatusCode DatabaseSetOption(struct AdbcDatabase* database, const char* key,
                                 const char* value, struct AdbcError* error) {
  SetError(error, "AdbcDatabaseSetOption not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode DatabaseSetOptionBytes(struct AdbcDatabase* database, const char* key,
                                      const uint8_t* value, size_t length,
                                      struct AdbcError* error) {
  SetError(error, "AdbcDatabaseSetOptionBytes not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode DatabaseSetOptionInt(struct AdbcDatabase* database, const char* key,
                                    int64_t value, struct AdbcError* error) {
  SetError(error, "AdbcDatabaseSetOptionInt not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode DatabaseSetOptionDouble(struct AdbcDatabase* database, const char* key,
                                       double value, struct AdbcError* error) {
  SetError(error, "AdbcDatabaseSetOptionDouble not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionCancel(struct AdbcConnection* connection,
                                struct AdbcError* error) {
  SetError(error, "AdbcConnectionCancel not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionCommit(struct AdbcConnection*, struct AdbcError* error) {
  SetError(error, "AdbcConnectionCommit not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionGetInfo(struct AdbcConnection* connection,
                                 const uint32_t* info_codes, size_t info_codes_length,
                                 struct ArrowArrayStream* out, struct AdbcError* error) {
  SetError(error, "AdbcConnectionGetInfo not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionGetObjects(struct AdbcConnection*, int, const char*, const char*,
                                    const char*, const char**, const char*,
                                    struct ArrowArrayStream*, struct AdbcError* error) {
  SetError(error, "AdbcConnectionGetObjects not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionGetOption(struct AdbcConnection* connection, const char* key,
                                   char* value, size_t* length, struct AdbcError* error) {
  SetError(error, "AdbcConnectionGetOption not implemented");
  return ADBC_STATUS_NOT_FOUND;
}

AdbcStatusCode ConnectionGetOptionBytes(struct AdbcConnection* connection,
                                        const char* key, uint8_t* value, size_t* length,
                                        struct AdbcError* error) {
  SetError(error, "AdbcConnectionGetOptionBytes not implemented");
  return ADBC_STATUS_NOT_FOUND;
}

AdbcStatusCode ConnectionGetOptionInt(struct AdbcConnection* connection, const char* key,
                                      int64_t* value, struct AdbcError* error) {
  SetError(error, "AdbcConnectionGetOptionInt not implemented");
  return ADBC_STATUS_NOT_FOUND;
}

AdbcStatusCode ConnectionGetOptionDouble(struct AdbcConnection* connection,
                                         const char* key, double* value,
                                         struct AdbcError* error) {
  SetError(error, "AdbcConnectionGetOptionDouble not implemented");
  return ADBC_STATUS_NOT_FOUND;
}

AdbcStatusCode ConnectionGetStatistics(struct AdbcConnection*, const char*, const char*,
                                       const char*, char, struct ArrowArrayStream*,
                                       struct AdbcError* error) {
  SetError(error, "AdbcConnectionGetStatistics not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionGetStatisticNames(struct AdbcConnection*,
                                           struct ArrowArrayStream*,
                                           struct AdbcError* error) {
  SetError(error, "AdbcConnectionGetStatisticNames not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionGetTableSchema(struct AdbcConnection*, const char*, const char*,
                                        const char*, struct ArrowSchema*,
                                        struct AdbcError* error) {
  SetError(error, "AdbcConnectionGetTableSchema not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionGetTableTypes(struct AdbcConnection*, struct ArrowArrayStream*,
                                       struct AdbcError* error) {
  SetError(error, "AdbcConnectionGetTableTypes not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionReadPartition(struct AdbcConnection* connection,
                                       const uint8_t* serialized_partition,
                                       size_t serialized_length,
                                       struct ArrowArrayStream* out,
                                       struct AdbcError* error) {
  SetError(error, "AdbcConnectionReadPartition not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionRollback(struct AdbcConnection*, struct AdbcError* error) {
  SetError(error, "AdbcConnectionRollback not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionSetOption(struct AdbcConnection*, const char*, const char*,
                                   struct AdbcError* error) {
  SetError(error, "AdbcConnectionSetOption not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionSetOptionBytes(struct AdbcConnection*, const char*,
                                        const uint8_t*, size_t, struct AdbcError* error) {
  SetError(error, "AdbcConnectionSetOptionBytes not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionSetOptionInt(struct AdbcConnection* connection, const char* key,
                                      int64_t value, struct AdbcError* error) {
  SetError(error, "AdbcConnectionSetOptionInt not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode ConnectionSetOptionDouble(struct AdbcConnection* connection,
                                         const char* key, double value,
                                         struct AdbcError* error) {
  SetError(error, "AdbcConnectionSetOptionDouble not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode StatementBind(struct AdbcStatement*, struct ArrowArray*,
                             struct ArrowSchema*, struct AdbcError* error) {
  SetError(error, "AdbcStatementBind not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode StatementBindStream(struct AdbcStatement*, struct ArrowArrayStream*,
                                   struct AdbcError* error) {
  SetError(error, "AdbcStatementBindStream not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode StatementCancel(struct AdbcStatement* statement, struct AdbcError* error) {
  SetError(error, "AdbcStatementCancel not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode StatementExecutePartitions(struct AdbcStatement* statement,
                                          struct ArrowSchema* schema,
                                          struct AdbcPartitions* partitions,
                                          int64_t* rows_affected,
                                          struct AdbcError* error) {
  SetError(error, "AdbcStatementExecutePartitions not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode StatementExecuteSchema(struct AdbcStatement* statement,
                                      struct ArrowSchema* schema,
                                      struct AdbcError* error) {
  SetError(error, "AdbcStatementExecuteSchema not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode StatementGetOption(struct AdbcStatement* statement, const char* key,
                                  char* value, size_t* length, struct AdbcError* error) {
  SetError(error, "AdbcStatementGetOption not implemented");
  return ADBC_STATUS_NOT_FOUND;
}

AdbcStatusCode StatementGetOptionBytes(struct AdbcStatement* statement, const char* key,
                                       uint8_t* value, size_t* length,
                                       struct AdbcError* error) {
  SetError(error, "AdbcStatementGetOptionBytes not implemented");
  return ADBC_STATUS_NOT_FOUND;
}

AdbcStatusCode StatementGetOptionInt(struct AdbcStatement* statement, const char* key,
                                     int64_t* value, struct AdbcError* error) {
  SetError(error, "AdbcStatementGetOptionInt not implemented");
  return ADBC_STATUS_NOT_FOUND;
}

AdbcStatusCode StatementGetOptionDouble(struct AdbcStatement* statement, const char* key,
                                        double* value, struct AdbcError* error) {
  SetError(error, "AdbcStatementGetOptionDouble not implemented");
  return ADBC_STATUS_NOT_FOUND;
}

AdbcStatusCode StatementGetParameterSchema(struct AdbcStatement* statement,
                                           struct ArrowSchema* schema,
                                           struct AdbcError* error) {
  SetError(error, "AdbcStatementGetParameterSchema not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode StatementPrepare(struct AdbcStatement*, struct AdbcError* error) {
  SetError(error, "AdbcStatementPrepare not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode StatementSetOption(struct AdbcStatement*, const char*, const char*,
                                  struct AdbcError* error) {
  SetError(error, "AdbcStatementSetOption not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode StatementSetOptionBytes(struct AdbcStatement*, const char*, const uint8_t*,
                                       size_t, struct AdbcError* error) {
  SetError(error, "AdbcStatementSetOptionBytes not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode StatementSetOptionInt(struct AdbcStatement* statement, const char* key,
                                     int64_t value, struct AdbcError* error) {
  SetError(error, "AdbcStatementSetOptionInt not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode StatementSetOptionDouble(struct AdbcStatement* statement, const char* key,
                                        double value, struct AdbcError* error) {
  SetError(error, "AdbcStatementSetOptionDouble not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode StatementSetSqlQuery(struct AdbcStatement*, const char*,
                                    struct AdbcError* error) {
  SetError(error, "AdbcStatementSetSqlQuery not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

AdbcStatusCode StatementSetSubstraitPlan(struct AdbcStatement*, const uint8_t*, size_t,
                                         struct AdbcError* error) {
  SetError(error, "AdbcStatementSetSubstraitPlan not implemented");
  return ADBC_STATUS_NOT_IMPLEMENTED;
}

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
  std::string additional_search_path_list;
};

/// Temporary state while the database is being configured.
struct TempConnection {
  std::unordered_map<std::string, std::string> options;
  std::unordered_map<std::string, std::string> bytes_options;
  std::unordered_map<std::string, int64_t> int_options;
  std::unordered_map<std::string, double> double_options;
};

static const char kDefaultEntrypoint[] = "AdbcDriverInit";
}  // namespace

// Other helpers (intentionally not in an anonymous namespace so they can be tested)
ADBC_EXPORT
std::filesystem::path InternalAdbcUserConfigDir() {
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
    config_dir /= "ADBC/Drivers";
  }
#elif defined(__APPLE__)
  auto dir = std::getenv("HOME");
  if (dir) {
    config_dir = std::filesystem::path(dir);
    config_dir /= "Library/Application Support/ADBC/Drivers";
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
    config_dir = config_dir / "adbc" / "drivers";
  }
#endif  // defined(_WIN32)

  return config_dir;
}

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
std::string InternalAdbcDriverManagerDefaultEntrypoint(const std::string& driver) {
  /// - libadbc_driver_sqlite.so.2.0.0 -> AdbcDriverSqliteInit
  /// - adbc_driver_sqlite.dll -> AdbcDriverSqliteInit
  /// - proprietary_driver.dll -> AdbcProprietaryDriverInit

  // Potential path -> filename
  // Treat both \ and / as directory separators on all platforms for simplicity
  std::string filename;
  {
    size_t pos = driver.find_last_of("/\\");
    if (pos != std::string::npos) {
      filename = driver.substr(pos + 1);
    } else {
      filename = driver;
    }
  }

  // Remove all extensions
  {
    size_t pos = filename.find('.');
    if (pos != std::string::npos) {
      filename = filename.substr(0, pos);
    }
  }

  // Remove lib prefix
  // https://stackoverflow.com/q/1878001/262727
  if (filename.rfind("lib", 0) == 0) {
    filename = filename.substr(3);
  }

  // Split on underscores, hyphens
  // Capitalize and join
  std::string entrypoint;
  entrypoint.reserve(filename.size());
  size_t pos = 0;
  while (pos < filename.size()) {
    size_t prev = pos;
    pos = filename.find_first_of("-_", pos);
    // if pos == npos this is the entire filename
    std::string token = filename.substr(prev, pos - prev);
    // capitalize first letter
    token[0] = std::toupper(static_cast<unsigned char>(token[0]));

    entrypoint += token;

    if (pos != std::string::npos) {
      pos++;
    }
  }

  if (entrypoint.rfind("Adbc", 0) != 0) {
    entrypoint = "Adbc" + entrypoint;
  }
  entrypoint += "Init";

  return entrypoint;
}

// Direct implementations of API methods

int AdbcErrorGetDetailCount(const struct AdbcError* error) {
  if (error->vendor_code == ADBC_ERROR_VENDOR_CODE_PRIVATE_DATA && error->private_data &&
      error->private_driver && error->private_driver->ErrorGetDetailCount) {
    return error->private_driver->ErrorGetDetailCount(error);
  }
  return 0;
}

struct AdbcErrorDetail AdbcErrorGetDetail(const struct AdbcError* error, int index) {
  if (error->vendor_code == ADBC_ERROR_VENDOR_CODE_PRIVATE_DATA && error->private_data &&
      error->private_driver && error->private_driver->ErrorGetDetail) {
    return error->private_driver->ErrorGetDetail(error, index);
  }
  return {nullptr, nullptr, 0};
}

const struct AdbcError* AdbcErrorFromArrayStream(struct ArrowArrayStream* stream,
                                                 AdbcStatusCode* status) {
  if (!stream->private_data || stream->release != ErrorArrayStreamRelease) {
    return nullptr;
  }
  auto* private_data = reinterpret_cast<struct ErrorArrayStream*>(stream->private_data);
  auto* error =
      private_data->private_driver->ErrorFromArrayStream(&private_data->stream, status);
  if (error) {
    const_cast<struct AdbcError*>(error)->private_driver = private_data->private_driver;
  }
  return error;
}

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

AdbcStatusCode AdbcDatabaseNew(struct AdbcDatabase* database, struct AdbcError* error) {
  // Allocate a temporary structure to store options pre-Init
  database->private_data = new TempDatabase();
  database->private_driver = nullptr;
  return ADBC_STATUS_OK;
}

AdbcStatusCode AdbcDatabaseGetOption(struct AdbcDatabase* database, const char* key,
                                     char* value, size_t* length,
                                     struct AdbcError* error) {
  if (database->private_driver) {
    INIT_ERROR(error, database);
    return database->private_driver->DatabaseGetOption(database, key, value, length,
                                                       error);
  }
  const auto* args = reinterpret_cast<const TempDatabase*>(database->private_data);
  const std::string* result = nullptr;
  if (std::strcmp(key, "driver") == 0) {
    result = &args->driver;
  } else if (std::strcmp(key, "entrypoint") == 0) {
    result = &args->entrypoint;
  } else {
    const auto it = args->options.find(key);
    if (it == args->options.end()) {
      SetError(error, std::string("Option not found: ") + key);
      return ADBC_STATUS_NOT_FOUND;
    }
    result = &it->second;
  }

  if (*length <= result->size() + 1) {
    // Enough space
    std::memcpy(value, result->c_str(), result->size() + 1);
  }
  *length = result->size() + 1;
  return ADBC_STATUS_OK;
}

AdbcStatusCode AdbcDatabaseGetOptionBytes(struct AdbcDatabase* database, const char* key,
                                          uint8_t* value, size_t* length,
                                          struct AdbcError* error) {
  if (database->private_driver) {
    INIT_ERROR(error, database);
    return database->private_driver->DatabaseGetOptionBytes(database, key, value, length,
                                                            error);
  }
  const auto* args = reinterpret_cast<const TempDatabase*>(database->private_data);
  const auto it = args->bytes_options.find(key);
  if (it == args->options.end()) {
    SetError(error, std::string("Option not found: ") + key);
    return ADBC_STATUS_NOT_FOUND;
  }
  const std::string& result = it->second;

  if (*length <= result.size()) {
    // Enough space
    std::memcpy(value, result.c_str(), result.size());
  }
  *length = result.size();
  return ADBC_STATUS_OK;
}

AdbcStatusCode AdbcDatabaseGetOptionInt(struct AdbcDatabase* database, const char* key,
                                        int64_t* value, struct AdbcError* error) {
  if (database->private_driver) {
    INIT_ERROR(error, database);
    return database->private_driver->DatabaseGetOptionInt(database, key, value, error);
  }
  const auto* args = reinterpret_cast<const TempDatabase*>(database->private_data);
  const auto it = args->int_options.find(key);
  if (it == args->int_options.end()) {
    SetError(error, std::string("Option not found: ") + key);
    return ADBC_STATUS_NOT_FOUND;
  }
  *value = it->second;
  return ADBC_STATUS_OK;
}

AdbcStatusCode AdbcDatabaseGetOptionDouble(struct AdbcDatabase* database, const char* key,
                                           double* value, struct AdbcError* error) {
  if (database->private_driver) {
    INIT_ERROR(error, database);
    return database->private_driver->DatabaseGetOptionDouble(database, key, value, error);
  }
  const auto* args = reinterpret_cast<const TempDatabase*>(database->private_data);
  const auto it = args->double_options.find(key);
  if (it == args->double_options.end()) {
    SetError(error, std::string("Option not found: ") + key);
    return ADBC_STATUS_NOT_FOUND;
  }
  *value = it->second;
  return ADBC_STATUS_OK;
}

AdbcStatusCode AdbcDatabaseSetOption(struct AdbcDatabase* database, const char* key,
                                     const char* value, struct AdbcError* error) {
  if (database->private_driver) {
    INIT_ERROR(error, database);
    return database->private_driver->DatabaseSetOption(database, key, value, error);
  }

  TempDatabase* args = reinterpret_cast<TempDatabase*>(database->private_data);
  if (std::strcmp(key, "driver") == 0) {
    args->driver = value;
  } else if (std::strcmp(key, "entrypoint") == 0) {
    args->entrypoint = value;
  } else {
    args->options[key] = value;
  }
  return ADBC_STATUS_OK;
}

AdbcStatusCode AdbcDatabaseSetOptionBytes(struct AdbcDatabase* database, const char* key,
                                          const uint8_t* value, size_t length,
                                          struct AdbcError* error) {
  if (database->private_driver) {
    INIT_ERROR(error, database);
    return database->private_driver->DatabaseSetOptionBytes(database, key, value, length,
                                                            error);
  }

  TempDatabase* args = reinterpret_cast<TempDatabase*>(database->private_data);
  args->bytes_options[key] = std::string(reinterpret_cast<const char*>(value), length);
  return ADBC_STATUS_OK;
}

AdbcStatusCode AdbcDatabaseSetOptionInt(struct AdbcDatabase* database, const char* key,
                                        int64_t value, struct AdbcError* error) {
  if (database->private_driver) {
    INIT_ERROR(error, database);
    return database->private_driver->DatabaseSetOptionInt(database, key, value, error);
  }

  TempDatabase* args = reinterpret_cast<TempDatabase*>(database->private_data);
  args->int_options[key] = value;
  return ADBC_STATUS_OK;
}

AdbcStatusCode AdbcDatabaseSetOptionDouble(struct AdbcDatabase* database, const char* key,
                                           double value, struct AdbcError* error) {
  if (database->private_driver) {
    INIT_ERROR(error, database);
    return database->private_driver->DatabaseSetOptionDouble(database, key, value, error);
  }

  TempDatabase* args = reinterpret_cast<TempDatabase*>(database->private_data);
  args->double_options[key] = value;
  return ADBC_STATUS_OK;
}

AdbcStatusCode AdbcDriverManagerDatabaseSetLoadFlags(struct AdbcDatabase* database,
                                                     AdbcLoadFlags flags,
                                                     struct AdbcError* error) {
  if (database->private_driver) {
    SetError(error, "Cannot SetLoadFlags after AdbcDatabaseInit");
    return ADBC_STATUS_INVALID_STATE;
  }

  TempDatabase* args = reinterpret_cast<TempDatabase*>(database->private_data);
  args->load_flags = flags;
  return ADBC_STATUS_OK;
}

AdbcStatusCode AdbcDriverManagerDatabaseSetAdditionalSearchPathList(
    struct AdbcDatabase* database, const char* path_list, struct AdbcError* error) {
  if (database->private_driver) {
    SetError(error, "Cannot SetAdditionalSearchPathList after AdbcDatabaseInit");
    return ADBC_STATUS_INVALID_STATE;
  }

  TempDatabase* args = reinterpret_cast<TempDatabase*>(database->private_data);
  if (path_list) {
    args->additional_search_path_list.assign(path_list);
  } else {
    args->additional_search_path_list.clear();
  }
  return ADBC_STATUS_OK;
}

AdbcStatusCode AdbcDriverManagerDatabaseSetInitFunc(struct AdbcDatabase* database,
                                                    AdbcDriverInitFunc init_func,
                                                    struct AdbcError* error) {
  if (database->private_driver) {
    SetError(error, "Cannot SetInitFunc after AdbcDatabaseInit");
    return ADBC_STATUS_INVALID_STATE;
  }

  TempDatabase* args = reinterpret_cast<TempDatabase*>(database->private_data);
  args->init_func = init_func;
  return ADBC_STATUS_OK;
}

AdbcStatusCode AdbcDatabaseInit(struct AdbcDatabase* database, struct AdbcError* error) {
  if (!database->private_data) {
    SetError(error, "Must call AdbcDatabaseNew before AdbcDatabaseInit");
    return ADBC_STATUS_INVALID_STATE;
  }
  TempDatabase* args = reinterpret_cast<TempDatabase*>(database->private_data);
  if (args->init_func) {
    // Do nothing
  } else if (args->driver.empty()) {
    SetError(error, "Must provide 'driver' parameter");
    return ADBC_STATUS_INVALID_ARGUMENT;
  }

  database->private_driver = new AdbcDriver;
  std::memset(database->private_driver, 0, sizeof(AdbcDriver));
  AdbcStatusCode status;
  // So we don't confuse a driver into thinking it's initialized already
  database->private_data = nullptr;
  if (args->init_func) {
    status = AdbcLoadDriverFromInitFunc(args->init_func, ADBC_VERSION_1_1_0,
                                        database->private_driver, error);
  } else if (!args->entrypoint.empty()) {
    status = AdbcFindLoadDriver(args->driver.c_str(), args->entrypoint.c_str(),
                                ADBC_VERSION_1_1_0, args->load_flags,
                                args->additional_search_path_list.data(),
                                database->private_driver, error);
  } else {
    status = AdbcFindLoadDriver(
        args->driver.c_str(), nullptr, ADBC_VERSION_1_1_0, args->load_flags,
        args->additional_search_path_list.data(), database->private_driver, error);
  }

  if (status != ADBC_STATUS_OK) {
    // Restore private_data so it will be released by AdbcDatabaseRelease
    database->private_data = args;
    if (database->private_driver->release) {
      database->private_driver->release(database->private_driver, error);
    }
    delete database->private_driver;
    database->private_driver = nullptr;
    return status;
  }

  // Errors that occur during AdbcDatabaseXXX() refer to the driver via
  // the private_driver member; however, after we return we have released
  // the driver and inspecting the error might segfault. Here, we scope
  // the driver-produced error to this function and make a copy if necessary.
  OwnedError driver_error;

  status = database->private_driver->DatabaseNew(database, &driver_error.error);
  if (status != ADBC_STATUS_OK) {
    if (database->private_driver->release) {
      SetError(error, &driver_error.error);
      database->private_driver->release(database->private_driver, nullptr);
    }
    delete database->private_driver;
    database->private_driver = nullptr;
    return status;
  }
  auto options = std::move(args->options);
  auto bytes_options = std::move(args->bytes_options);
  auto int_options = std::move(args->int_options);
  auto double_options = std::move(args->double_options);
  delete args;

  INIT_ERROR(error, database);
  for (const auto& option : options) {
    status = database->private_driver->DatabaseSetOption(
        database, option.first.c_str(), option.second.c_str(), &driver_error.error);
    if (status != ADBC_STATUS_OK) break;
  }
  for (const auto& option : bytes_options) {
    status = database->private_driver->DatabaseSetOptionBytes(
        database, option.first.c_str(),
        reinterpret_cast<const uint8_t*>(option.second.data()), option.second.size(),
        &driver_error.error);
    if (status != ADBC_STATUS_OK) break;
  }
  for (const auto& option : int_options) {
    status = database->private_driver->DatabaseSetOptionInt(
        database, option.first.c_str(), option.second, &driver_error.error);
    if (status != ADBC_STATUS_OK) break;
  }
  for (const auto& option : double_options) {
    status = database->private_driver->DatabaseSetOptionDouble(
        database, option.first.c_str(), option.second, &driver_error.error);
    if (status != ADBC_STATUS_OK) break;
  }

  if (status != ADBC_STATUS_OK) {
    // Release the database
    std::ignore = database->private_driver->DatabaseRelease(database, nullptr);
    if (database->private_driver->release) {
      SetError(error, &driver_error.error);
      database->private_driver->release(database->private_driver, nullptr);
    }
    delete database->private_driver;
    database->private_driver = nullptr;
    // Should be redundant, but ensure that AdbcDatabaseRelease
    // below doesn't think that it contains a TempDatabase
    database->private_data = nullptr;
    return status;
  }

  return database->private_driver->DatabaseInit(database, error);
}

AdbcStatusCode AdbcDatabaseRelease(struct AdbcDatabase* database,
                                   struct AdbcError* error) {
  if (!database->private_driver) {
    if (database->private_data) {
      TempDatabase* args = reinterpret_cast<TempDatabase*>(database->private_data);
      delete args;
      database->private_data = nullptr;
      return ADBC_STATUS_OK;
    }
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, database);
  auto status = database->private_driver->DatabaseRelease(database, error);
  if (database->private_driver->release) {
    database->private_driver->release(database->private_driver, error);
  }
  delete database->private_driver;
  database->private_data = nullptr;
  database->private_driver = nullptr;
  return status;
}

AdbcStatusCode AdbcConnectionCancel(struct AdbcConnection* connection,
                                    struct AdbcError* error) {
  if (!connection->private_driver) {
    SetError(error, "AdbcConnectionCancel: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, connection);
  return connection->private_driver->ConnectionCancel(connection, error);
}

AdbcStatusCode AdbcConnectionCommit(struct AdbcConnection* connection,
                                    struct AdbcError* error) {
  if (!connection->private_driver) {
    SetError(error, "AdbcConnectionCommit: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, connection);
  return connection->private_driver->ConnectionCommit(connection, error);
}

AdbcStatusCode AdbcConnectionGetInfo(struct AdbcConnection* connection,
                                     const uint32_t* info_codes, size_t info_codes_length,
                                     struct ArrowArrayStream* out,
                                     struct AdbcError* error) {
  if (!connection->private_driver) {
    SetError(error, "AdbcConnectionGetInfo: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, connection);
  WRAP_STREAM(connection->private_driver->ConnectionGetInfo(
                  connection, info_codes, info_codes_length, out, error),
              out, connection);
}

AdbcStatusCode AdbcConnectionGetObjects(struct AdbcConnection* connection, int depth,
                                        const char* catalog, const char* db_schema,
                                        const char* table_name, const char** table_types,
                                        const char* column_name,
                                        struct ArrowArrayStream* stream,
                                        struct AdbcError* error) {
  if (!connection->private_driver) {
    SetError(error, "AdbcConnectionGetObjects: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, connection);
  WRAP_STREAM(connection->private_driver->ConnectionGetObjects(
                  connection, depth, catalog, db_schema, table_name, table_types,
                  column_name, stream, error),
              stream, connection);
}

AdbcStatusCode AdbcConnectionGetOption(struct AdbcConnection* connection, const char* key,
                                       char* value, size_t* length,
                                       struct AdbcError* error) {
  if (!connection->private_data) {
    SetError(error, "AdbcConnectionGetOption: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  if (!connection->private_driver) {
    // Init not yet called, get the saved option
    const auto* args = reinterpret_cast<const TempConnection*>(connection->private_data);
    const auto it = args->options.find(key);
    if (it == args->options.end()) {
      return ADBC_STATUS_NOT_FOUND;
    }
    if (*length >= it->second.size() + 1) {
      std::memcpy(value, it->second.c_str(), it->second.size() + 1);
    }
    *length = it->second.size() + 1;
    return ADBC_STATUS_OK;
  }
  INIT_ERROR(error, connection);
  return connection->private_driver->ConnectionGetOption(connection, key, value, length,
                                                         error);
}

AdbcStatusCode AdbcConnectionGetOptionBytes(struct AdbcConnection* connection,
                                            const char* key, uint8_t* value,
                                            size_t* length, struct AdbcError* error) {
  if (!connection->private_data) {
    SetError(error, "AdbcConnectionGetOption: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  if (!connection->private_driver) {
    // Init not yet called, get the saved option
    const auto* args = reinterpret_cast<const TempConnection*>(connection->private_data);
    const auto it = args->bytes_options.find(key);
    if (it == args->options.end()) {
      return ADBC_STATUS_NOT_FOUND;
    }
    if (*length >= it->second.size() + 1) {
      std::memcpy(value, it->second.data(), it->second.size() + 1);
    }
    *length = it->second.size() + 1;
    return ADBC_STATUS_OK;
  }
  INIT_ERROR(error, connection);
  return connection->private_driver->ConnectionGetOptionBytes(connection, key, value,
                                                              length, error);
}

AdbcStatusCode AdbcConnectionGetOptionInt(struct AdbcConnection* connection,
                                          const char* key, int64_t* value,
                                          struct AdbcError* error) {
  if (!connection->private_data) {
    SetError(error, "AdbcConnectionGetOption: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  if (!connection->private_driver) {
    // Init not yet called, get the saved option
    const auto* args = reinterpret_cast<const TempConnection*>(connection->private_data);
    const auto it = args->int_options.find(key);
    if (it == args->int_options.end()) {
      return ADBC_STATUS_NOT_FOUND;
    }
    *value = it->second;
    return ADBC_STATUS_OK;
  }
  INIT_ERROR(error, connection);
  return connection->private_driver->ConnectionGetOptionInt(connection, key, value,
                                                            error);
}

AdbcStatusCode AdbcConnectionGetOptionDouble(struct AdbcConnection* connection,
                                             const char* key, double* value,
                                             struct AdbcError* error) {
  if (!connection->private_data) {
    SetError(error, "AdbcConnectionGetOption: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  if (!connection->private_driver) {
    // Init not yet called, get the saved option
    const auto* args = reinterpret_cast<const TempConnection*>(connection->private_data);
    const auto it = args->double_options.find(key);
    if (it == args->double_options.end()) {
      return ADBC_STATUS_NOT_FOUND;
    }
    *value = it->second;
    return ADBC_STATUS_OK;
  }
  INIT_ERROR(error, connection);
  return connection->private_driver->ConnectionGetOptionDouble(connection, key, value,
                                                               error);
}

AdbcStatusCode AdbcConnectionGetStatistics(struct AdbcConnection* connection,
                                           const char* catalog, const char* db_schema,
                                           const char* table_name, char approximate,
                                           struct ArrowArrayStream* out,
                                           struct AdbcError* error) {
  if (!connection->private_driver) {
    SetError(error, "AdbcConnectionGetStatistics: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, connection);
  WRAP_STREAM(
      connection->private_driver->ConnectionGetStatistics(
          connection, catalog, db_schema, table_name, approximate == 1, out, error),
      out, connection);
}

AdbcStatusCode AdbcConnectionGetStatisticNames(struct AdbcConnection* connection,
                                               struct ArrowArrayStream* out,
                                               struct AdbcError* error) {
  if (!connection->private_driver) {
    SetError(error, "AdbcConnectionGetStatisticNames: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, connection);
  WRAP_STREAM(
      connection->private_driver->ConnectionGetStatisticNames(connection, out, error),
      out, connection);
}

AdbcStatusCode AdbcConnectionGetTableSchema(struct AdbcConnection* connection,
                                            const char* catalog, const char* db_schema,
                                            const char* table_name,
                                            struct ArrowSchema* schema,
                                            struct AdbcError* error) {
  if (!connection->private_driver) {
    SetError(error, "AdbcConnectionGetTableSchema: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, connection);
  return connection->private_driver->ConnectionGetTableSchema(
      connection, catalog, db_schema, table_name, schema, error);
}

AdbcStatusCode AdbcConnectionGetTableTypes(struct AdbcConnection* connection,
                                           struct ArrowArrayStream* stream,
                                           struct AdbcError* error) {
  if (!connection->private_driver) {
    SetError(error, "AdbcConnectionGetTableTypes: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, connection);
  WRAP_STREAM(
      connection->private_driver->ConnectionGetTableTypes(connection, stream, error),
      stream, connection);
}

AdbcStatusCode AdbcConnectionInit(struct AdbcConnection* connection,
                                  struct AdbcDatabase* database,
                                  struct AdbcError* error) {
  if (!connection->private_data) {
    SetError(error, "Must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  } else if (!database->private_driver) {
    SetError(error, "Database is not initialized");
    return ADBC_STATUS_INVALID_ARGUMENT;
  }
  TempConnection* args = reinterpret_cast<TempConnection*>(connection->private_data);
  connection->private_data = nullptr;
  std::unordered_map<std::string, std::string> options = std::move(args->options);
  std::unordered_map<std::string, std::string> bytes_options =
      std::move(args->bytes_options);
  std::unordered_map<std::string, int64_t> int_options = std::move(args->int_options);
  std::unordered_map<std::string, double> double_options =
      std::move(args->double_options);
  delete args;

  auto status = database->private_driver->ConnectionNew(connection, error);
  if (status != ADBC_STATUS_OK) return status;
  connection->private_driver = database->private_driver;

  for (const auto& option : options) {
    status = database->private_driver->ConnectionSetOption(
        connection, option.first.c_str(), option.second.c_str(), error);
    if (status != ADBC_STATUS_OK) return status;
  }
  for (const auto& option : bytes_options) {
    status = database->private_driver->ConnectionSetOptionBytes(
        connection, option.first.c_str(),
        reinterpret_cast<const uint8_t*>(option.second.data()), option.second.size(),
        error);
    if (status != ADBC_STATUS_OK) return status;
  }
  for (const auto& option : int_options) {
    status = database->private_driver->ConnectionSetOptionInt(
        connection, option.first.c_str(), option.second, error);
    if (status != ADBC_STATUS_OK) return status;
  }
  for (const auto& option : double_options) {
    status = database->private_driver->ConnectionSetOptionDouble(
        connection, option.first.c_str(), option.second, error);
    if (status != ADBC_STATUS_OK) return status;
  }
  INIT_ERROR(error, connection);
  return connection->private_driver->ConnectionInit(connection, database, error);
}

AdbcStatusCode AdbcConnectionNew(struct AdbcConnection* connection,
                                 struct AdbcError* error) {
  // Allocate a temporary structure to store options pre-Init, because
  // we don't get access to the database (and hence the driver
  // function table) until then
  connection->private_data = new TempConnection;
  connection->private_driver = nullptr;
  return ADBC_STATUS_OK;
}

AdbcStatusCode AdbcConnectionReadPartition(struct AdbcConnection* connection,
                                           const uint8_t* serialized_partition,
                                           size_t serialized_length,
                                           struct ArrowArrayStream* out,
                                           struct AdbcError* error) {
  if (!connection->private_driver) {
    SetError(error, "AdbcConnectionReadPartition: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, connection);
  WRAP_STREAM(connection->private_driver->ConnectionReadPartition(
                  connection, serialized_partition, serialized_length, out, error),
              out, connection);
}

AdbcStatusCode AdbcConnectionRelease(struct AdbcConnection* connection,
                                     struct AdbcError* error) {
  if (!connection->private_driver) {
    if (connection->private_data) {
      TempConnection* args = reinterpret_cast<TempConnection*>(connection->private_data);
      delete args;
      connection->private_data = nullptr;
      return ADBC_STATUS_OK;
    }
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, connection);
  auto status = connection->private_driver->ConnectionRelease(connection, error);
  connection->private_driver = nullptr;
  return status;
}

AdbcStatusCode AdbcConnectionRollback(struct AdbcConnection* connection,
                                      struct AdbcError* error) {
  if (!connection->private_driver) {
    SetError(error, "AdbcConnectionRollback: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, connection);
  return connection->private_driver->ConnectionRollback(connection, error);
}

AdbcStatusCode AdbcConnectionSetOption(struct AdbcConnection* connection, const char* key,
                                       const char* value, struct AdbcError* error) {
  if (!connection->private_data) {
    SetError(error, "AdbcConnectionSetOption: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  if (!connection->private_driver) {
    // Init not yet called, save the option
    TempConnection* args = reinterpret_cast<TempConnection*>(connection->private_data);
    args->options[key] = value;
    return ADBC_STATUS_OK;
  }
  INIT_ERROR(error, connection);
  return connection->private_driver->ConnectionSetOption(connection, key, value, error);
}

AdbcStatusCode AdbcConnectionSetOptionBytes(struct AdbcConnection* connection,
                                            const char* key, const uint8_t* value,
                                            size_t length, struct AdbcError* error) {
  if (!connection->private_data) {
    SetError(error, "AdbcConnectionSetOptionInt: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  if (!connection->private_driver) {
    // Init not yet called, save the option
    TempConnection* args = reinterpret_cast<TempConnection*>(connection->private_data);
    args->bytes_options[key] = std::string(reinterpret_cast<const char*>(value), length);
    return ADBC_STATUS_OK;
  }
  INIT_ERROR(error, connection);
  return connection->private_driver->ConnectionSetOptionBytes(connection, key, value,
                                                              length, error);
}

AdbcStatusCode AdbcConnectionSetOptionInt(struct AdbcConnection* connection,
                                          const char* key, int64_t value,
                                          struct AdbcError* error) {
  if (!connection->private_data) {
    SetError(error, "AdbcConnectionSetOptionInt: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  if (!connection->private_driver) {
    // Init not yet called, save the option
    TempConnection* args = reinterpret_cast<TempConnection*>(connection->private_data);
    args->int_options[key] = value;
    return ADBC_STATUS_OK;
  }
  INIT_ERROR(error, connection);
  return connection->private_driver->ConnectionSetOptionInt(connection, key, value,
                                                            error);
}

AdbcStatusCode AdbcConnectionSetOptionDouble(struct AdbcConnection* connection,
                                             const char* key, double value,
                                             struct AdbcError* error) {
  if (!connection->private_data) {
    SetError(error, "AdbcConnectionSetOptionDouble: must call AdbcConnectionNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  if (!connection->private_driver) {
    // Init not yet called, save the option
    TempConnection* args = reinterpret_cast<TempConnection*>(connection->private_data);
    args->double_options[key] = value;
    return ADBC_STATUS_OK;
  }
  INIT_ERROR(error, connection);
  return connection->private_driver->ConnectionSetOptionDouble(connection, key, value,
                                                               error);
}

AdbcStatusCode AdbcStatementBind(struct AdbcStatement* statement,
                                 struct ArrowArray* values, struct ArrowSchema* schema,
                                 struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementBind: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementBind(statement, values, schema, error);
}

AdbcStatusCode AdbcStatementBindStream(struct AdbcStatement* statement,
                                       struct ArrowArrayStream* stream,
                                       struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementBindStream: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementBindStream(statement, stream, error);
}

AdbcStatusCode AdbcStatementCancel(struct AdbcStatement* statement,
                                   struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementCancel: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementCancel(statement, error);
}

// XXX: cpplint gets confused here if declared as 'struct ArrowSchema* schema'
AdbcStatusCode AdbcStatementExecutePartitions(struct AdbcStatement* statement,
                                              ArrowSchema* schema,
                                              struct AdbcPartitions* partitions,
                                              int64_t* rows_affected,
                                              struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementExecutePartitions: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementExecutePartitions(
      statement, schema, partitions, rows_affected, error);
}

AdbcStatusCode AdbcStatementExecuteQuery(struct AdbcStatement* statement,
                                         struct ArrowArrayStream* out,
                                         int64_t* rows_affected,
                                         struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementExecuteQuery: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  WRAP_STREAM(statement->private_driver->StatementExecuteQuery(statement, out,
                                                               rows_affected, error),
              out, statement);
}

AdbcStatusCode AdbcStatementExecuteSchema(struct AdbcStatement* statement,
                                          struct ArrowSchema* schema,
                                          struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementExecuteSchema: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementExecuteSchema(statement, schema, error);
}

AdbcStatusCode AdbcStatementGetOption(struct AdbcStatement* statement, const char* key,
                                      char* value, size_t* length,
                                      struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementGetOption: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementGetOption(statement, key, value, length,
                                                       error);
}

AdbcStatusCode AdbcStatementGetOptionBytes(struct AdbcStatement* statement,
                                           const char* key, uint8_t* value,
                                           size_t* length, struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementGetOptionBytes: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementGetOptionBytes(statement, key, value, length,
                                                            error);
}

AdbcStatusCode AdbcStatementGetOptionInt(struct AdbcStatement* statement, const char* key,
                                         int64_t* value, struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementGetOptionInt: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementGetOptionInt(statement, key, value, error);
}

AdbcStatusCode AdbcStatementGetOptionDouble(struct AdbcStatement* statement,
                                            const char* key, double* value,
                                            struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementGetOptionDouble: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementGetOptionDouble(statement, key, value,
                                                             error);
}

AdbcStatusCode AdbcStatementGetParameterSchema(struct AdbcStatement* statement,
                                               struct ArrowSchema* schema,
                                               struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementGetParameterSchema: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementGetParameterSchema(statement, schema, error);
}

AdbcStatusCode AdbcStatementNew(struct AdbcConnection* connection,
                                struct AdbcStatement* statement,
                                struct AdbcError* error) {
  if (!connection->private_driver) {
    SetError(error, "AdbcStatementNew: must call AdbcConnectionInit first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, connection);
  auto status = connection->private_driver->StatementNew(connection, statement, error);
  statement->private_driver = connection->private_driver;
  return status;
}

AdbcStatusCode AdbcStatementPrepare(struct AdbcStatement* statement,
                                    struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementPrepare: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementPrepare(statement, error);
}

AdbcStatusCode AdbcStatementRelease(struct AdbcStatement* statement,
                                    struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementRelease: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  auto status = statement->private_driver->StatementRelease(statement, error);
  statement->private_driver = nullptr;
  return status;
}

AdbcStatusCode AdbcStatementSetOption(struct AdbcStatement* statement, const char* key,
                                      const char* value, struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementSetOption: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementSetOption(statement, key, value, error);
}

AdbcStatusCode AdbcStatementSetOptionBytes(struct AdbcStatement* statement,
                                           const char* key, const uint8_t* value,
                                           size_t length, struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementSetOptionBytes: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementSetOptionBytes(statement, key, value, length,
                                                            error);
}

AdbcStatusCode AdbcStatementSetOptionInt(struct AdbcStatement* statement, const char* key,
                                         int64_t value, struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementSetOptionInt: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementSetOptionInt(statement, key, value, error);
}

AdbcStatusCode AdbcStatementSetOptionDouble(struct AdbcStatement* statement,
                                            const char* key, double value,
                                            struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementSetOptionDouble: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementSetOptionDouble(statement, key, value,
                                                             error);
}

AdbcStatusCode AdbcStatementSetSqlQuery(struct AdbcStatement* statement,
                                        const char* query, struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementSetSqlQuery: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementSetSqlQuery(statement, query, error);
}

AdbcStatusCode AdbcStatementSetSubstraitPlan(struct AdbcStatement* statement,
                                             const uint8_t* plan, size_t length,
                                             struct AdbcError* error) {
  if (!statement->private_driver) {
    SetError(error, "AdbcStatementSetSubstraitPlan: must call AdbcStatementNew first");
    return ADBC_STATUS_INVALID_STATE;
  }
  INIT_ERROR(error, statement);
  return statement->private_driver->StatementSetSubstraitPlan(statement, plan, length,
                                                              error);
}

const char* AdbcStatusCodeMessage(AdbcStatusCode code) {
#define CASE(CONSTANT)         \
  case ADBC_STATUS_##CONSTANT: \
    return #CONSTANT;

  switch (code) {
    CASE(OK);
    CASE(UNKNOWN);
    CASE(NOT_IMPLEMENTED);
    CASE(NOT_FOUND);
    CASE(ALREADY_EXISTS);
    CASE(INVALID_ARGUMENT);
    CASE(INVALID_STATE);
    CASE(INVALID_DATA);
    CASE(INTEGRITY);
    CASE(INTERNAL);
    CASE(IO);
    CASE(CANCELLED);
    CASE(TIMEOUT);
    CASE(UNAUTHENTICATED);
    CASE(UNAUTHORIZED);
    default:
      return "(invalid code)";
  }
#undef CASE
}

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

AdbcStatusCode AdbcLoadDriver(const char* driver_name, const char* entrypoint,
                              int version, void* raw_driver, struct AdbcError* error) {
  // maintain old behavior of allowing relative paths (because dlopen allows it)
  // but don't enable searching for manifests by default. It will need to be explicitly
  // enabled by calling AdbcFindLoadDriver directly.
  return AdbcFindLoadDriver(driver_name, entrypoint, version,
                            ADBC_LOAD_FLAG_ALLOW_RELATIVE_PATHS, nullptr, raw_driver,
                            error);
}

AdbcStatusCode AdbcLoadDriverFromInitFunc(AdbcDriverInitFunc init_func, int version,
                                          void* raw_driver, struct AdbcError* error) {
  constexpr std::array<int, 2> kSupportedVersions = {
      ADBC_VERSION_1_1_0,
      ADBC_VERSION_1_0_0,
  };

  if (!raw_driver) {
    SetError(error, "Must provide non-NULL raw_driver");
    return ADBC_STATUS_INVALID_ARGUMENT;
  }

  switch (version) {
    case ADBC_VERSION_1_0_0:
    case ADBC_VERSION_1_1_0:
      break;
    default:
      SetError(error, "Only ADBC 1.0.0 and 1.1.0 are supported");
      return ADBC_STATUS_NOT_IMPLEMENTED;
  }

#define FILL_DEFAULT(DRIVER, STUB) \
  if (!DRIVER->STUB) {             \
    DRIVER->STUB = &STUB;          \
  }
#define CHECK_REQUIRED(DRIVER, STUB)                                           \
  if (!DRIVER->STUB) {                                                         \
    SetError(error, "Driver does not implement required function Adbc" #STUB); \
    return ADBC_STATUS_INTERNAL;                                               \
  }

  // Starting from the passed version, try each (older) version in
  // succession with the underlying driver until we find one that's
  // accepted.
  AdbcStatusCode result = ADBC_STATUS_NOT_IMPLEMENTED;
  for (const int try_version : kSupportedVersions) {
    if (try_version > version) continue;
    result = init_func(try_version, raw_driver, error);
    if (result != ADBC_STATUS_NOT_IMPLEMENTED) break;
  }
  if (result != ADBC_STATUS_OK) {
    return result;
  }

  if (version >= ADBC_VERSION_1_0_0) {
    auto* driver = reinterpret_cast<struct AdbcDriver*>(raw_driver);
    CHECK_REQUIRED(driver, DatabaseNew);
    CHECK_REQUIRED(driver, DatabaseInit);
    CHECK_REQUIRED(driver, DatabaseRelease);
    FILL_DEFAULT(driver, DatabaseSetOption);

    CHECK_REQUIRED(driver, ConnectionNew);
    CHECK_REQUIRED(driver, ConnectionInit);
    CHECK_REQUIRED(driver, ConnectionRelease);
    FILL_DEFAULT(driver, ConnectionCommit);
    FILL_DEFAULT(driver, ConnectionGetInfo);
    FILL_DEFAULT(driver, ConnectionGetObjects);
    FILL_DEFAULT(driver, ConnectionGetTableSchema);
    FILL_DEFAULT(driver, ConnectionGetTableTypes);
    FILL_DEFAULT(driver, ConnectionReadPartition);
    FILL_DEFAULT(driver, ConnectionRollback);
    FILL_DEFAULT(driver, ConnectionSetOption);

    FILL_DEFAULT(driver, StatementExecutePartitions);
    CHECK_REQUIRED(driver, StatementExecuteQuery);
    CHECK_REQUIRED(driver, StatementNew);
    CHECK_REQUIRED(driver, StatementRelease);
    FILL_DEFAULT(driver, StatementBind);
    FILL_DEFAULT(driver, StatementBindStream);
    FILL_DEFAULT(driver, StatementGetParameterSchema);
    FILL_DEFAULT(driver, StatementPrepare);
    FILL_DEFAULT(driver, StatementSetOption);
    FILL_DEFAULT(driver, StatementSetSqlQuery);
    FILL_DEFAULT(driver, StatementSetSubstraitPlan);
  }
  if (version >= ADBC_VERSION_1_1_0) {
    auto* driver = reinterpret_cast<struct AdbcDriver*>(raw_driver);
    FILL_DEFAULT(driver, ErrorGetDetailCount);
    FILL_DEFAULT(driver, ErrorGetDetail);
    FILL_DEFAULT(driver, ErrorFromArrayStream);

    FILL_DEFAULT(driver, DatabaseGetOption);
    FILL_DEFAULT(driver, DatabaseGetOptionBytes);
    FILL_DEFAULT(driver, DatabaseGetOptionDouble);
    FILL_DEFAULT(driver, DatabaseGetOptionInt);
    FILL_DEFAULT(driver, DatabaseSetOptionBytes);
    FILL_DEFAULT(driver, DatabaseSetOptionDouble);
    FILL_DEFAULT(driver, DatabaseSetOptionInt);

    FILL_DEFAULT(driver, ConnectionCancel);
    FILL_DEFAULT(driver, ConnectionGetOption);
    FILL_DEFAULT(driver, ConnectionGetOptionBytes);
    FILL_DEFAULT(driver, ConnectionGetOptionDouble);
    FILL_DEFAULT(driver, ConnectionGetOptionInt);
    FILL_DEFAULT(driver, ConnectionGetStatistics);
    FILL_DEFAULT(driver, ConnectionGetStatisticNames);
    FILL_DEFAULT(driver, ConnectionSetOptionBytes);
    FILL_DEFAULT(driver, ConnectionSetOptionDouble);
    FILL_DEFAULT(driver, ConnectionSetOptionInt);

    FILL_DEFAULT(driver, StatementCancel);
    FILL_DEFAULT(driver, StatementExecuteSchema);
    FILL_DEFAULT(driver, StatementGetOption);
    FILL_DEFAULT(driver, StatementGetOptionBytes);
    FILL_DEFAULT(driver, StatementGetOptionDouble);
    FILL_DEFAULT(driver, StatementGetOptionInt);
    FILL_DEFAULT(driver, StatementSetOptionBytes);
    FILL_DEFAULT(driver, StatementSetOptionDouble);
    FILL_DEFAULT(driver, StatementSetOptionInt);
  }

  return ADBC_STATUS_OK;

#undef FILL_DEFAULT
#undef CHECK_REQUIRED
}
