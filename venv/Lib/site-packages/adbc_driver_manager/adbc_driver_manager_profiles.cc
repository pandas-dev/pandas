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
#endif                // defined(_WIN32)

#include "adbc_driver_manager_internal.h"

#include <filesystem>
#include <regex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <toml++/toml.hpp>

#include "arrow-adbc/adbc.h"
#include "arrow-adbc/adbc_driver_manager.h"

using namespace std::string_literals;  // NOLINT [build/namespaces]

namespace {

#ifdef _WIN32
static const wchar_t* kAdbcProfilePath = L"ADBC_PROFILE_PATH";
#else
static const char* kAdbcProfilePath = "ADBC_PROFILE_PATH";
#endif  // _WIN32

}  // namespace

// FilesystemProfile needs external linkage for use in internal header
struct FilesystemProfile {
  std::filesystem::path path;
  std::string driver;
  std::unordered_map<std::string, std::string> options;
  std::unordered_map<std::string, int64_t> int_options;
  std::unordered_map<std::string, double> double_options;

  std::vector<const char*> options_keys;
  std::vector<const char*> options_values;

  std::vector<const char*> int_option_keys;
  std::vector<int64_t> int_option_values;

  std::vector<const char*> double_option_keys;
  std::vector<double> double_option_values;

  void PopulateConnectionProfile(struct AdbcConnectionProfile* out) {
    options_keys.reserve(options.size());
    options_values.reserve(options.size());
    for (const auto& [key, value] : options) {
      options_keys.push_back(key.c_str());
      options_values.push_back(value.c_str());
    }

    int_option_keys.reserve(int_options.size());
    int_option_values.reserve(int_options.size());
    for (const auto& [key, value] : int_options) {
      int_option_keys.push_back(key.c_str());
      int_option_values.push_back(value);
    }

    double_option_keys.reserve(double_options.size());
    double_option_values.reserve(double_options.size());
    for (const auto& [key, value] : double_options) {
      double_option_keys.push_back(key.c_str());
      double_option_values.push_back(value);
    }

    out->private_data = new FilesystemProfile(std::move(*this));
    out->release = [](AdbcConnectionProfile* profile) {
      if (!profile || !profile->private_data) {
        return;
      }

      delete static_cast<FilesystemProfile*>(profile->private_data);
      profile->private_data = nullptr;
      profile->release = nullptr;
    };

    out->GetDriverName = [](AdbcConnectionProfile* profile, const char** out,
                            AdbcDriverInitFunc* init_func,
                            struct AdbcError* error) -> AdbcStatusCode {
      if (!profile || !profile->private_data) {
        SetError(error, "Invalid connection profile");
        return ADBC_STATUS_INVALID_ARGUMENT;
      }

      auto* fs_profile = static_cast<FilesystemProfile*>(profile->private_data);
      *out = fs_profile->driver.c_str();
      *init_func = nullptr;
      return ADBC_STATUS_OK;
    };

    out->GetOptions = [](AdbcConnectionProfile* profile, const char*** keys,
                         const char*** values, size_t* num_options,
                         struct AdbcError* error) -> AdbcStatusCode {
      if (!profile || !profile->private_data) {
        SetError(error, "Invalid connection profile");
        return ADBC_STATUS_INVALID_ARGUMENT;
      }

      if (!keys || !values || !num_options) {
        SetError(error, "Output parameters cannot be null");
        return ADBC_STATUS_INVALID_ARGUMENT;
      }

      auto* fs_profile = static_cast<FilesystemProfile*>(profile->private_data);
      *num_options = fs_profile->options.size();
      *keys = fs_profile->options_keys.data();
      *values = fs_profile->options_values.data();
      return ADBC_STATUS_OK;
    };

    out->GetIntOptions = [](AdbcConnectionProfile* profile, const char*** keys,
                            const int64_t** values, size_t* num_options,
                            struct AdbcError* error) -> AdbcStatusCode {
      if (!profile || !profile->private_data) {
        SetError(error, "Invalid connection profile");
        return ADBC_STATUS_INVALID_ARGUMENT;
      }

      if (!keys || !values || !num_options) {
        SetError(error, "Output parameters cannot be null");
        return ADBC_STATUS_INVALID_ARGUMENT;
      }

      auto* fs_profile = static_cast<FilesystemProfile*>(profile->private_data);
      *num_options = fs_profile->int_options.size();
      *keys = fs_profile->int_option_keys.data();
      *values = fs_profile->int_option_values.data();
      return ADBC_STATUS_OK;
    };

    out->GetDoubleOptions = [](AdbcConnectionProfile* profile, const char*** keys,
                               const double** values, size_t* num_options,
                               struct AdbcError* error) -> AdbcStatusCode {
      if (!profile || !profile->private_data) {
        SetError(error, "Invalid connection profile");
        return ADBC_STATUS_INVALID_ARGUMENT;
      }

      if (!keys || !values || !num_options) {
        SetError(error, "Output parameters cannot be null");
        return ADBC_STATUS_INVALID_ARGUMENT;
      }

      auto* fs_profile = static_cast<FilesystemProfile*>(profile->private_data);
      *num_options = fs_profile->double_options.size();
      *keys = fs_profile->double_option_keys.data();
      *values = fs_profile->double_option_values.data();
      return ADBC_STATUS_OK;
    };
  }
};

struct ProfileVisitor {
  FilesystemProfile& profile;
  const std::filesystem::path& profile_path;
  struct AdbcError* error;

  bool VisitTable(const std::string& prefix, toml::table& table) {
    for (const auto& [key, value] : table) {
      if (auto* str = value.as_string()) {
        profile.options[prefix + key.data()] = str->get();
      } else if (auto* int_val = value.as_integer()) {
        profile.int_options[prefix + key.data()] = int_val->get();
      } else if (auto* double_val = value.as_floating_point()) {
        profile.double_options[prefix + key.data()] = double_val->get();
      } else if (auto* bool_val = value.as_boolean()) {
        profile.options[prefix + key.data()] = bool_val->get() ? "true" : "false";
      } else if (value.is_table()) {
        if (!VisitTable(prefix + key.data() + ".", *value.as_table())) {
          return false;
        }
      } else {
        std::string message = "Unsupported value type for key '" +
                              std::string(key.str()) + "' in profile '" +
                              profile_path.string() + "'";
        SetError(error, std::move(message));
        return false;
      }
    }
    return !error->message;
  }
};

// Public implementations (non-static for use across translation units)
AdbcStatusCode ProcessProfileValue(std::string_view key, std::string_view value,
                                   std::string& out, struct AdbcError* error) {
  if (value.empty()) {
    out = "";
    return ADBC_STATUS_OK;
  }

  static const std::regex pattern(R"(\{\{\s*([^{}]*?)\s*\}\})");
  auto end_of_last_match = value.begin();
  auto begin = std::regex_iterator(value.begin(), value.end(), pattern);
  auto end = decltype(begin){};
  std::match_results<std::string_view::iterator>::difference_type pos_last_match = 0;

  out.resize(0);
  for (auto itr = begin; itr != end; ++itr) {
    auto match = *itr;
    auto pos_match = match.position();
    auto diff = pos_match - pos_last_match;
    auto start_match = end_of_last_match;
    std::advance(start_match, diff);
    out.append(end_of_last_match, start_match);

    const auto content = match[1].str();
    if (content.rfind("env_var(", 0) != 0) {
      std::string message = "In profile: unsupported interpolation type in key `" +
                            std::string(key) + "`: `" + content + "`";
      SetError(error, message);
      return ADBC_STATUS_INVALID_ARGUMENT;
    }

    if (content[content.size() - 1] != ')') {
      std::string message = "In profile: malformed env_var() in key `" +
                            std::string(key) + "`: missing closing parenthesis";
      SetError(error, message);
      return ADBC_STATUS_INVALID_ARGUMENT;
    }

    const auto env_var_name = content.substr(8, content.size() - 9);
    if (env_var_name.empty()) {
      std::string message = "In profile: malformed env_var() in key `" +
                            std::string(key) + "`: missing environment variable name";
      SetError(error, message);
      return ADBC_STATUS_INVALID_ARGUMENT;
    }

    std::string env_var_value;
#ifdef _WIN32
    auto local_env_var = Utf8Decode(std::string(env_var_name));
    DWORD required_size = GetEnvironmentVariableW(local_env_var.c_str(), NULL, 0);
    if (required_size != 0) {
      std::wstring wvalue;
      wvalue.resize(required_size);
      DWORD actual_size =
          GetEnvironmentVariableW(local_env_var.c_str(), wvalue.data(), required_size);
      // remove null terminator
      wvalue.resize(actual_size);
      env_var_value = Utf8Encode(wvalue);
    }
#else
    const char* env_value = std::getenv(env_var_name.c_str());
    if (env_value) {
      env_var_value = std::string(env_value);
    }
#endif
    out.append(env_var_value);

    auto length_match = match.length();
    pos_last_match = pos_match + length_match;
    end_of_last_match = start_match;
    std::advance(end_of_last_match, length_match);
  }

  out.append(end_of_last_match, value.end());
  return ADBC_STATUS_OK;
}

AdbcStatusCode LoadProfileFile(const std::filesystem::path& profile_path,
                               FilesystemProfile& profile, struct AdbcError* error) {
  toml::table config;
  try {
    config = toml::parse_file(profile_path.native());
  } catch (const toml::parse_error& err) {
    std::string message = "Could not open profile. ";
    message += err.what();
    message += ". Profile: ";
    message += profile_path.string();
    SetError(error, std::move(message));
    return ADBC_STATUS_INVALID_ARGUMENT;
  }

  profile.path = profile_path;
  if (!config["profile_version"].is_integer()) {
    std::string message =
        "Profile version is not an integer in profile '" + profile_path.string() + "'";
    SetError(error, std::move(message));
    return ADBC_STATUS_INVALID_ARGUMENT;
  }

  const auto version = config["profile_version"].value_or(int64_t(1));
  switch (version) {
    case 1:
      break;
    default: {
      std::string message =
          "Profile version '" + std::to_string(version) +
          "' is not supported by this driver manager. Profile: " + profile_path.string();
      SetError(error, std::move(message));
      return ADBC_STATUS_INVALID_ARGUMENT;
    }
  }

  profile.driver = config["driver"].value_or(""s);

  auto options = config.at_path("Options");
  if (!options.is_table()) {
    std::string message =
        "Profile options is not a table in profile '" + profile_path.string() + "'";
    SetError(error, std::move(message));
    return ADBC_STATUS_INVALID_ARGUMENT;
  }

  auto* options_table = options.as_table();
  ProfileVisitor v{profile, profile_path, error};
  if (!v.VisitTable("", *options_table)) {
    return ADBC_STATUS_INVALID_ARGUMENT;
  }

  return ADBC_STATUS_OK;
}

namespace {

static SearchPaths GetProfileSearchPaths(const char* additional_search_path_list) {
  SearchPaths search_paths;
  {
    std::vector<std::filesystem::path> additional_paths;
    if (additional_search_path_list) {
      additional_paths = InternalAdbcParsePath(additional_search_path_list);
    }

    for (const auto& path : additional_paths) {
      search_paths.emplace_back(SearchPathSource::kAdditional, path);
    }
  }

  {
    auto env_paths = GetEnvPaths(kAdbcProfilePath);
    search_paths.insert(search_paths.end(), env_paths.begin(), env_paths.end());
  }

#if ADBC_CONDA_BUILD
#ifdef _WIN32
  const wchar_t* conda_name = L"CONDA_PREFIX";
#else
  const char* conda_name = "CONDA_PREFIX";
#endif  // _WIN32

  auto venv = GetEnvPaths(conda_name);
  for (const auto& [_, venv_path] : venv) {
    search_paths.emplace_back(SearchPathSource::kConda,
                              venv_path / "etc" / "adbc" / "profiles");
  }
#else
  search_paths.emplace_back(SearchPathSource::kDisabledAtCompileTime, "Conda prefix");
#endif  // ADBC_CONDA_BUILD

#ifdef _WIN32
  const wchar_t* profiles_dir = L"Profiles";
#elif defined(__APPLE__)
  const char* profiles_dir = "Profiles";
#else
  const char* profiles_dir = "profiles";
#endif  // defined(_WIN32)

  auto user_dir = InternalAdbcUserConfigDir() / profiles_dir;
  search_paths.emplace_back(SearchPathSource::kUser, user_dir);
  return search_paths;
}

}  // namespace

AdbcStatusCode AdbcProfileProviderFilesystem(const char* profile_name,
                                             const char* additional_search_path_list,
                                             struct AdbcConnectionProfile* out,
                                             struct AdbcError* error) {
  if (profile_name == nullptr || strlen(profile_name) == 0) {
    SetError(error, "Profile name is empty");
    return ADBC_STATUS_INVALID_ARGUMENT;
  }

  if (!out) {
    SetError(error, "Output profile is null");
    return ADBC_STATUS_INVALID_ARGUMENT;
  }

  std::memset(out, 0, sizeof(*out));
  std::filesystem::path profile_path(profile_name);
  if (profile_path.has_extension()) {
    if (HasExtension(profile_path, ".toml")) {
      if (!std::filesystem::exists(profile_path)) {
        SetError(error, "Profile file does not exist: " + profile_path.string());
        return ADBC_STATUS_NOT_FOUND;
      }

      FilesystemProfile profile;
      CHECK_STATUS(LoadProfileFile(profile_path, profile, error));
      profile.PopulateConnectionProfile(out);
      return ADBC_STATUS_OK;
    }
  }

  if (profile_path.is_absolute()) {
    profile_path.replace_extension(".toml");

    FilesystemProfile profile;
    CHECK_STATUS(LoadProfileFile(profile_path, profile, error));
    profile.PopulateConnectionProfile(out);
    return ADBC_STATUS_OK;
  }

  SearchPaths search_paths = GetProfileSearchPaths(additional_search_path_list);
  SearchPaths extra_debug_info;
  for (const auto& [source, search_path] : search_paths) {
    if (source == SearchPathSource::kRegistry || source == SearchPathSource::kUnset ||
        source == SearchPathSource::kDoesNotExist ||
        source == SearchPathSource::kDisabledAtCompileTime ||
        source == SearchPathSource::kDisabledAtRunTime ||
        source == SearchPathSource::kOtherError) {
      continue;
    }

    std::filesystem::path full_path = search_path / profile_path;
    full_path.replace_extension(".toml");
    if (std::filesystem::exists(full_path)) {
      OwnedError intermediate_error;

      FilesystemProfile profile;
      auto status = LoadProfileFile(full_path, profile, &intermediate_error.error);
      if (status == ADBC_STATUS_OK) {
        profile.PopulateConnectionProfile(out);
        return ADBC_STATUS_OK;
      } else if (status == ADBC_STATUS_INVALID_ARGUMENT) {
        search_paths.insert(search_paths.end(), extra_debug_info.begin(),
                            extra_debug_info.end());
        if (intermediate_error.error.message) {
          std::string error_message = intermediate_error.error.message;
          // Remove [Driver Manager] prefix so it doesn't get repeated
          error_message = error_message.substr(17);
          AddSearchPathsToError(search_paths, SearchPathType::kProfile, error_message);
          SetError(error, error_message);
        }
        return status;
      }

      std::string message = "found ";
      message += full_path.string();
      message += " but: ";
      if (intermediate_error.error.message) {
        std::string m = intermediate_error.error.message;
        // Remove [Driver Manager] prefix so it doesn't get repeated
        message += m.substr(17);
      } else {
        message += "could not load the profile";
      }
      extra_debug_info.emplace_back(SearchPathSource::kOtherError, std::move(message));
    }
  }

  search_paths.insert(search_paths.end(), extra_debug_info.begin(),
                      extra_debug_info.end());
  std::string error_message = "Profile not found: " + std::string(profile_name);
  AddSearchPathsToError(search_paths, SearchPathType::kProfile, error_message);
  SetError(error, std::move(error_message));
  return ADBC_STATUS_NOT_FOUND;
}
