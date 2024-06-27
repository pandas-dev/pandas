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

#include "adbc_driver_manager.h"
#include <adbc.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <string>
#include <unordered_map>
#include <utility>

#if defined(_WIN32)
#include <windows.h>  // Must come first

#include <libloaderapi.h>
#include <strsafe.h>
#else
#include <dlfcn.h>
#endif  // defined(_WIN32)

namespace {

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

// Driver state

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
    // TODO(apache/arrow-adbc#204): causes tests to segfault
    // Need to refcount the driver DLL; also, errors may retain a reference to
    // release() from the DLL - how to handle this?
  }

  AdbcStatusCode Load(const char* library, struct AdbcError* error) {
    std::string error_message;
#if defined(_WIN32)
    HMODULE handle = LoadLibraryExA(library, NULL, 0);
    if (!handle) {
      error_message += library;
      error_message += ": LoadLibraryExA() failed: ";
      GetWinError(&error_message);

      std::string full_driver_name = library;
      full_driver_name += ".dll";
      handle = LoadLibraryExA(full_driver_name.c_str(), NULL, 0);
      if (!handle) {
        error_message += '\n';
        error_message += full_driver_name;
        error_message += ": LoadLibraryExA() failed: ";
        GetWinError(&error_message);
      }
    }
    if (!handle) {
      SetError(error, error_message);
      return ADBC_STATUS_INTERNAL;
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
      error_message = "[DriverManager] dlopen() failed: ";
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
      SetError(error, error_message);
      return ADBC_STATUS_INTERNAL;
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
std::string AdbcDriverManagerDefaultEntrypoint(const std::string& driver) {
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
      error->private_driver) {
    return error->private_driver->ErrorGetDetailCount(error);
  }
  return 0;
}

struct AdbcErrorDetail AdbcErrorGetDetail(const struct AdbcError* error, int index) {
  if (error->vendor_code == ADBC_ERROR_VENDOR_CODE_PRIVATE_DATA && error->private_data &&
      error->private_driver) {
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
    status = AdbcLoadDriver(args->driver.c_str(), args->entrypoint.c_str(),
                            ADBC_VERSION_1_1_0, database->private_driver, error);
  } else {
    status = AdbcLoadDriver(args->driver.c_str(), nullptr, ADBC_VERSION_1_1_0,
                            database->private_driver, error);
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
  status = database->private_driver->DatabaseNew(database, error);
  if (status != ADBC_STATUS_OK) {
    if (database->private_driver->release) {
      database->private_driver->release(database->private_driver, error);
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
    status = database->private_driver->DatabaseSetOption(database, option.first.c_str(),
                                                         option.second.c_str(), error);
    if (status != ADBC_STATUS_OK) break;
  }
  for (const auto& option : bytes_options) {
    status = database->private_driver->DatabaseSetOptionBytes(
        database, option.first.c_str(),
        reinterpret_cast<const uint8_t*>(option.second.data()), option.second.size(),
        error);
    if (status != ADBC_STATUS_OK) break;
  }
  for (const auto& option : int_options) {
    status = database->private_driver->DatabaseSetOptionInt(
        database, option.first.c_str(), option.second, error);
    if (status != ADBC_STATUS_OK) break;
  }
  for (const auto& option : double_options) {
    status = database->private_driver->DatabaseSetOptionDouble(
        database, option.first.c_str(), option.second, error);
    if (status != ADBC_STATUS_OK) break;
  }

  if (status != ADBC_STATUS_OK) {
    // Release the database
    std::ignore = database->private_driver->DatabaseRelease(database, error);
    if (database->private_driver->release) {
      database->private_driver->release(database->private_driver, error);
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

AdbcStatusCode AdbcLoadDriver(const char* driver_name, const char* entrypoint,
                              int version, void* raw_driver, struct AdbcError* error) {
  AdbcDriverInitFunc init_func;
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
    SetError(error, "Must provide non-NULL raw_driver");
    return ADBC_STATUS_INVALID_ARGUMENT;
  }
  auto* driver = reinterpret_cast<struct AdbcDriver*>(raw_driver);

  ManagedLibrary library;
  AdbcStatusCode status = library.Load(driver_name, error);
  if (status != ADBC_STATUS_OK) {
    // AdbcDatabaseInit tries to call this if set
    driver->release = nullptr;
    return status;
  }

  void* load_handle = nullptr;
  if (entrypoint) {
    status = library.Lookup(entrypoint, &load_handle, error);
  } else {
    auto name = AdbcDriverManagerDefaultEntrypoint(driver_name);
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
