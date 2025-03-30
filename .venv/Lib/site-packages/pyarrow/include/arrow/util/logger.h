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

#include <chrono>
#include <iosfwd>
#include <memory>
#include <string_view>

#include "arrow/result.h"
#include "arrow/status.h"
#include "arrow/util/logging.h"
#include "arrow/util/macros.h"
#include "arrow/util/visibility.h"

namespace arrow {
namespace util {

struct SourceLocation {
  const char* file = "";
  int line = 0;
};

struct LogDetails {
  ArrowLogLevel severity = ArrowLogLevel::ARROW_INFO;
  std::chrono::system_clock::time_point timestamp = std::chrono::system_clock::now();
  SourceLocation source_location{};
  std::string_view message = "";
};

/// \brief A base interface for custom loggers.
///
/// Loggers can be added to the LoggerRegistry for global access or directly provided to
/// certain logging utilities.
class Logger {
 public:
  virtual ~Logger() = default;

  virtual void Log(const LogDetails& details) = 0;

  virtual bool Flush(std::chrono::microseconds timeout) { return true; }
  bool Flush() { return this->Flush(std::chrono::microseconds::max()); }

  virtual bool is_enabled() const { return true; }

  virtual ArrowLogLevel severity_threshold() const { return ArrowLogLevel::ARROW_TRACE; }
};

/// \brief Creates a simple logger that redirects output to std::cerr
ARROW_EXPORT std::shared_ptr<Logger> MakeOStreamLogger(ArrowLogLevel severity_threshold);
/// \brief Creates a simple logger that redirects output to the provided ostream
ARROW_EXPORT std::shared_ptr<Logger> MakeOStreamLogger(ArrowLogLevel severity_threshold,
                                                       std::ostream& sink);

class ARROW_EXPORT LoggerRegistry {
 public:
  /// \brief Add a logger to the registry with the associated name
  ///
  /// Returns Invalid if a logger with the provided name already exists. Users should call
  /// `UnregisterLogger` first if they wish to overwrite it.
  static Status RegisterLogger(std::string_view name, std::shared_ptr<Logger> logger);

  /// \brief Remove a logger from the registry
  static void UnregisterLogger(std::string_view name);

  /// \brief Return the logger associated with the provided name
  ///
  /// If `name` is empty, the default logger is returned. If `name` doesn't match any of
  /// the registered loggers then a non-null noop logger is returned
  static std::shared_ptr<Logger> GetLogger(std::string_view name = "");

  /// \brief Return the default logger
  static std::shared_ptr<Logger> GetDefaultLogger();
  /// \brief Set the default logger
  static void SetDefaultLogger(std::shared_ptr<Logger> logger);
};

/// \brief Represents a single log record to be emitted by an underlying logger
class ARROW_EXPORT LogMessage {
 public:
  /// \brief Construct a LogMessage with the provided underlying logger
  LogMessage(ArrowLogLevel severity, std::shared_ptr<Logger> logger,
             SourceLocation source_location = {});
  /// \brief Construct a LogMessage with the provided logger name, which will be used to
  /// find an underlying logger in the registry
  LogMessage(ArrowLogLevel severity, std::string_view logger_name,
             SourceLocation source_location = {});

  std::ostream& Stream();

  // Convenience method - mainly for use in ARROW_LOG_* macros. This prevents unnecessary
  // argument evaluation when log statements are stripped in certain builds
  template <typename... Args>
  LogMessage& Append(Args&&... args) {
    if constexpr (sizeof...(Args) > 0) {
      if (CheckIsEnabled()) {
        (Stream() << ... << args);
      }
    }
    return *this;
  }

 private:
  bool CheckIsEnabled();

  class Impl;
  std::shared_ptr<Impl> impl_;
};

}  // namespace util
}  // namespace arrow

// For the following macros, log statements with a lower severity than
// `ARROW_MINIMUM_LOG_LEVEL` will be stripped from the build
#ifndef ARROW_MINIMUM_LOG_LEVEL
#  define ARROW_MINIMUM_LOG_LEVEL -1000
#endif

#define ARROW_LOGGER_INTERNAL(LOGGER, LEVEL)                                      \
  (::arrow::util::LogMessage(::arrow::util::ArrowLogLevel::ARROW_##LEVEL, LOGGER, \
                             ::arrow::util::SourceLocation{__FILE__, __LINE__}))

static_assert(static_cast<int>(::arrow::util::ArrowLogLevel::ARROW_TRACE) == -2);
#if ARROW_MINIMUM_LOG_LEVEL <= -2
#  define ARROW_LOGGER_TRACE(LOGGER, ...) \
    (ARROW_LOGGER_INTERNAL(LOGGER, TRACE).Append(__VA_ARGS__))
#else
#  define ARROW_LOGGER_TRACE(...) ARROW_UNUSED(0)
#endif

static_assert(static_cast<int>(::arrow::util::ArrowLogLevel::ARROW_DEBUG) == -1);
#if ARROW_MINIMUM_LOG_LEVEL <= -1
#  define ARROW_LOGGER_DEBUG(LOGGER, ...) \
    (ARROW_LOGGER_INTERNAL(LOGGER, DEBUG).Append(__VA_ARGS__))
#else
#  define ARROW_LOGGER_DEBUG(...) ARROW_UNUSED(0)
#endif

static_assert(static_cast<int>(::arrow::util::ArrowLogLevel::ARROW_INFO) == 0);
#if ARROW_MINIMUM_LOG_LEVEL <= 0
#  define ARROW_LOGGER_INFO(LOGGER, ...) \
    (ARROW_LOGGER_INTERNAL(LOGGER, INFO).Append(__VA_ARGS__))
#else
#  define ARROW_LOGGER_INFO(...) ARROW_UNUSED(0)
#endif

static_assert(static_cast<int>(::arrow::util::ArrowLogLevel::ARROW_WARNING) == 1);
#if ARROW_MINIMUM_LOG_LEVEL <= 1
#  define ARROW_LOGGER_WARNING(LOGGER, ...) \
    (ARROW_LOGGER_INTERNAL(LOGGER, WARNING).Append(__VA_ARGS__))
#else
#  define ARROW_LOGGER_WARNING(...) ARROW_UNUSED(0)
#endif

static_assert(static_cast<int>(::arrow::util::ArrowLogLevel::ARROW_ERROR) == 2);
#if ARROW_MINIMUM_LOG_LEVEL <= 2
#  define ARROW_LOGGER_ERROR(LOGGER, ...) \
    (ARROW_LOGGER_INTERNAL(LOGGER, ERROR).Append(__VA_ARGS__))
#else
#  define ARROW_LOGGER_ERROR(...) ARROW_UNUSED(0)
#endif

static_assert(static_cast<int>(::arrow::util::ArrowLogLevel::ARROW_FATAL) == 3);
#if ARROW_MINIMUM_LOG_LEVEL <= 3
#  define ARROW_LOGGER_FATAL(LOGGER, ...) \
    (ARROW_LOGGER_INTERNAL(LOGGER, FATAL).Append(__VA_ARGS__))
#else
#  define ARROW_LOGGER_FATAL(...) ARROW_UNUSED(0)
#endif

#define ARROW_LOGGER_CALL(LOGGER, LEVEL, ...) ARROW_LOGGER_##LEVEL(LOGGER, __VA_ARGS__)
