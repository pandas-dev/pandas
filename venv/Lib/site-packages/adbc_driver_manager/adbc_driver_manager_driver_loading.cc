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

// Driver loading and resolution implementation for ADBC driver manager

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>  // Must come first

#ifndef NTDDI_VERSION
#define NTDDI_VERSION 0x0A00000C  // For SHGetKnownFolderPath in ShlObj_core.h in ShlObj.h
#endif

#include <libloaderapi.h>
#include <string.h>  // _wcsnicmp
#else
#include <dlfcn.h>
#endif  // defined(_WIN32)

#include <toml++/toml.hpp>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "adbc_driver_manager_internal.h"
#include "arrow-adbc/adbc.h"
#include "arrow-adbc/adbc_driver_manager.h"
#include "current_arch.h"

using namespace std::string_literals;  // NOLINT [build/namespaces]

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

/// \brief Parses a TOML manifest file and extracts driver information.
///
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
    }
    std::string message = "Driver path not found in manifest '";
    message += driver_manifest.string();
    message += "' for current architecture '";
    message += adbc::CurrentArch();
    message += "'. Value was not a string";
    SetError(error, std::move(message));
    return ADBC_STATUS_INVALID_ARGUMENT;
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
  DWORD required_size = GetEnvironmentVariableW(env_var, NULL, 0);
  if (required_size == 0) {
    return {};
  }

  std::wstring path_var;
  path_var.resize(required_size);
  DWORD actual_size = GetEnvironmentVariableW(env_var, path_var.data(), required_size);
  // Remove null terminator
  path_var.resize(actual_size);
  auto path = Utf8Encode(path_var);
#else
  const char* path_var = std::getenv(env_var);
  if (!path_var) {
    return {};
  }
  std::string path(path_var);
#endif  // _WIN32
  SearchPaths paths;
  for (auto parsed_path : InternalAdbcParsePath(path)) {
    paths.emplace_back(SearchPathSource::kEnv, parsed_path);
  }
  return paths;
}

#if defined(_WIN32) || defined(__APPLE__)
static constexpr const char* kDriversSubdir = "Drivers";
#else
static constexpr const char* kDriversSubdir = "drivers";
#endif

SearchPaths GetSearchPaths(const AdbcLoadFlags levels) {
  SearchPaths paths;
  if (levels & ADBC_LOAD_FLAG_SEARCH_ENV) {
    // Check the ADBC_DRIVER_PATH environment variable
    paths = GetEnvPaths(kAdbcDriverPath);
  }

  if (levels & ADBC_LOAD_FLAG_SEARCH_USER) {
    // Check the user configuration directory
    std::filesystem::path user_config_dir = InternalAdbcUserConfigDir() / kDriversSubdir;
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
#if !defined(_WIN32)
    const std::filesystem::path system_config_dir = InternalAdbcSystemConfigDir();
    if (std::filesystem::exists(system_config_dir)) {
      paths.emplace_back(SearchPathSource::kSystem, std::move(system_config_dir));
    } else {
      paths.emplace_back(SearchPathSource::kDoesNotExist, std::move(system_config_dir));
    }
#endif  // !defined(_WIN32)
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

void ManagedLibrary::Release() {
  // TODO(apache/arrow-adbc#204): causes tests to segfault.  Need to
  // refcount the driver DLL; also, errors may retain a reference to
  // release() from the DLL - how to handle this?  It's unlikely we can
  // actually do this - in general shared libraries are not safe to unload.
}

/// \brief Resolve the driver name to a concrete location.
AdbcStatusCode ManagedLibrary::GetDriverInfo(
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
        return Load(info.lib_path.native(), {}, error);
      }
      return status;
    }

    // if the extension is not .toml, then just try to load the provided
    // path as if it was an absolute path to a driver library
    info.lib_path = driver_path;
    return Load(driver_path.native(), {}, error);
  }

  if (driver_path.is_absolute()) {
    // if we have an absolute path without an extension, first see if there's a
    // toml file with the same name.
    driver_path.replace_extension(".toml");
    if (std::filesystem::exists(driver_path)) {
      auto status = LoadDriverManifest(driver_path, info, error);
      if (status == ADBC_STATUS_OK) {
        return Load(info.lib_path.native(), {}, error);
      }
    }

    driver_path.replace_extension("");
    // otherwise just try to load the provided path as if it was an absolute path
    info.lib_path = driver_path;
    return Load(driver_path.native(), {}, error);
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
      return Load(driver_path.native(), {}, error);
    }

    SetError(error, "Driver name has unrecognized extension: " +
                        driver_path.extension().string());
    return ADBC_STATUS_INVALID_ARGUMENT;
  }

  // not an absolute path, no extension. Let's search the configured paths
  // based on the options
  // FindDriver will set info.lib_path
  // XXX(lidavidm): the control flow in this call chain is excessively
  // convoluted and it's hard to determine if DriverInfo is fully
  // initialized or not in all non-error paths
  return FindDriver(driver_path, load_options, additional_search_paths, info, error);
}

/// \return ADBC_STATUS_NOT_FOUND if the driver shared library could not be
///   found (via dlopen) or if a manifest was found but did not contain a
///   path for the current platform, ADBC_STATUS_INVALID_ARGUMENT if a
///   manifest was found but could not be parsed, ADBC_STATUS_OK otherwise
///
/// May modify search_paths to add error info
AdbcStatusCode ManagedLibrary::SearchPathsForDriver(
    const std::filesystem::path& driver_path, SearchPaths& search_paths, DriverInfo& info,
    struct AdbcError* error) {
  SearchPaths extra_debug_info;
  for (const auto& [source, search_path] : search_paths) {
    if (source == SearchPathSource::kRegistry || source == SearchPathSource::kUnset ||
        source == SearchPathSource::kDoesNotExist ||
        source == SearchPathSource::kDisabledAtCompileTime ||
        source == SearchPathSource::kDisabledAtRunTime ||
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
        status = Load(info.lib_path.native(), {}, &intermediate_error.error);
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
        extra_debug_info.emplace_back(SearchPathSource::kOtherError, std::move(message));
        search_paths.insert(search_paths.end(), extra_debug_info.begin(),
                            extra_debug_info.end());
        return status;
      } else if (status == ADBC_STATUS_INVALID_ARGUMENT) {
        // The manifest was invalid. Don't ignore that!
        search_paths.insert(search_paths.end(), extra_debug_info.begin(),
                            extra_debug_info.end());
        if (intermediate_error.error.message) {
          std::string error_message = intermediate_error.error.message;
          AddSearchPathsToError(search_paths, SearchPathType::kManifest, error_message);
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
    auto status = Load(full_path.native(), {}, nullptr);
    if (status == ADBC_STATUS_OK) {
      info.lib_path = full_path;
      return status;
    }
  }

  search_paths.insert(search_paths.end(), extra_debug_info.begin(),
                      extra_debug_info.end());
  return ADBC_STATUS_NOT_FOUND;
}

AdbcStatusCode ManagedLibrary::FindDriver(
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
    if (!(load_options & ADBC_LOAD_FLAG_SEARCH_ENV)) {
      search_paths.emplace_back(SearchPathSource::kDisabledAtRunTime,
                                "ADBC_DRIVER_PATH (enable ADBC_LOAD_FLAG_SEARCH_ENV)");
    } else if (search_paths.empty()) {
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
        for (const auto& [_, venv_path] : venv) {
          search_paths.emplace_back(SearchPathSource::kConda,
                                    venv_path / "etc" / "adbc" / "drivers");
        }
      }
    } else {
      search_paths.emplace_back(SearchPathSource::kDisabledAtRunTime,
                                "Conda prefix (enable ADBC_LOAD_FLAG_SEARCH_ENV)");
    }
#else
    if (load_options & ADBC_LOAD_FLAG_SEARCH_ENV) {
      search_paths.emplace_back(SearchPathSource::kDisabledAtCompileTime, "Conda prefix");
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
      return Load(info.lib_path.native(), {}, error);
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
  } else {
    search_paths.emplace_back(SearchPathSource::kDisabledAtRunTime,
                              "HKEY_CURRENT_USER (enable ADBC_LOAD_FLAG_SEARCH_USER)");
  }

  if (load_options & ADBC_LOAD_FLAG_SEARCH_SYSTEM) {
    // Check the system registry for the driver.
    auto status =
        LoadDriverFromRegistry(HKEY_LOCAL_MACHINE, driver_path.native(), info, error);
    if (status == ADBC_STATUS_OK) {
      return Load(info.lib_path.native(), {}, error);
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
  } else {
    search_paths.emplace_back(SearchPathSource::kDisabledAtRunTime,
                              "HKEY_LOCAL_MACHINE (enable ADBC_LOAD_FLAG_SEARCH_SYSTEM)");
  }

  info.lib_path = driver_path;
  return Load(driver_path.native(), search_paths, error);
#else
  // Otherwise, search the configured paths.
  SearchPaths more_search_paths =
      GetSearchPaths(load_options & ~ADBC_LOAD_FLAG_SEARCH_ENV);
  auto status = SearchPathsForDriver(driver_path, more_search_paths, info, error);
  if (status == ADBC_STATUS_NOT_FOUND) {
    if (!(load_options & ADBC_LOAD_FLAG_SEARCH_USER)) {
      std::filesystem::path user_config_dir =
          InternalAdbcUserConfigDir() / kDriversSubdir;
      std::string message = "user config dir ";
      message += user_config_dir.string();
      message += " (enable ADBC_LOAD_FLAG_SEARCH_USER)";
      more_search_paths.emplace_back(SearchPathSource::kDisabledAtRunTime,
                                     std::move(message));
    }
    // Windows searches registry keys, so this only applies to other OSes
#if !defined(_WIN32)
    if (!(load_options & ADBC_LOAD_FLAG_SEARCH_SYSTEM)) {
      std::filesystem::path system_config_dir = InternalAdbcSystemConfigDir();
      std::string message = "system config dir ";
      message += system_config_dir.string();
      message += " (enable ADBC_LOAD_FLAG_SEARCH_SYSTEM)";
      more_search_paths.emplace_back(SearchPathSource::kDisabledAtRunTime,
                                     std::move(message));
    }
#endif  // !defined(_WIN32)

    // If we reach here, we didn't find the driver in any of the paths
    // so let's just attempt to load it as default behavior.
    search_paths.insert(search_paths.end(), more_search_paths.begin(),
                        more_search_paths.end());
    info.lib_path = driver_path;
    return Load(driver_path.native(), search_paths, error);
  }
  return status;
#endif  // _WIN32
}

/// \return ADBC_STATUS_NOT_FOUND if the driver shared library could not be
///   found, ADBC_STATUS_OK otherwise
AdbcStatusCode ManagedLibrary::Load(const string_type& library,
                                    const SearchPaths& attempted_paths,
                                    struct AdbcError* error) {
  std::string error_message;
#if defined(_WIN32)
  HMODULE handle = LoadLibraryExW(library.c_str(), NULL, 0);
  if (!handle) {
    error_message = "Could not load `";
    error_message += Utf8Encode(library);
    error_message += "`: LoadLibraryExW() failed: ";
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
    std::string name = Utf8Encode(library);
    std::string message = CheckNonPrintableLibraryName(name);
    if (!message.empty()) {
      error_message += "\n";
      error_message += message;
    }
    AddSearchPathsToError(attempted_paths, SearchPathType::kManifest, error_message);
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

  void* handle = dlopen(library.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (!handle) {
    error_message = "Could not load `";
    error_message += library;
    error_message += "`: dlopen() failed: ";
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
                           kPlatformLibrarySuffix.size(), kPlatformLibrarySuffix) != 0) {
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
    std::string message = CheckNonPrintableLibraryName(library);
    if (!message.empty()) {
      error_message += "\n";
      error_message += message;
    }
    AddSearchPathsToError(attempted_paths, SearchPathType::kManifest, error_message);
    SetError(error, error_message);
    return ADBC_STATUS_NOT_FOUND;
  }
#endif  // defined(_WIN32)
  return ADBC_STATUS_OK;
}

AdbcStatusCode ManagedLibrary::Lookup(const char* name, void** func,
                                      struct AdbcError* error) {
#if defined(_WIN32)
  void* load_handle = reinterpret_cast<void*>(GetProcAddress(handle, name));
  if (!load_handle) {
    std::string message = "GetProcAddress(";
    message += name;
    message += ") failed: ";
    GetWinError(&message);
    AppendError(error, message);
    return ADBC_STATUS_INTERNAL;
  }
#else
  void* load_handle = dlsym(handle, name);
  if (!load_handle) {
    std::string message = "dlsym(";
    message += name;
    message += ") failed: ";
    message += dlerror();
    AppendError(error, message);
    return ADBC_STATUS_INTERNAL;
  }
#endif  // defined(_WIN32)
  *func = load_handle;
  return ADBC_STATUS_OK;
}

ADBC_EXPORT
std::string InternalAdbcDriverManagerDefaultEntrypoint(const std::string& driver) {
  /// - libadbc_driver_sqlite.so.2.0.0 -> AdbcDriverSqliteInit
  /// - adbc_driver_sqlite.dll -> AdbcDriverSqliteInit
  /// - proprietary_driver.dll -> AdbcProprietaryDriverInit

  // N.B.(https://github.com/apache/arrow-adbc/issues/3680): sanity checks
  assert(!driver.empty());

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
    token[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(token[0])));

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
