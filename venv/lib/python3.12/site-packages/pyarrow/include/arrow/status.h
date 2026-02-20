// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.
//
// A Status encapsulates the result of an operation.  It may indicate success,
// or it may indicate an error with an associated error message.
//
// Multiple threads can invoke const methods on a Status without
// external synchronization, but if any of the threads may call a
// non-const method, all threads accessing the same Status must use
// external synchronization.

// Adapted from Apache Kudu, TensorFlow

#pragma once

#include <cstring>
#include <iosfwd>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

#include "arrow/util/compare.h"
#include "arrow/util/macros.h"
#include "arrow/util/string_util.h"
#include "arrow/util/visibility.h"

#ifdef ARROW_EXTRA_ERROR_CONTEXT

/// \brief Return with given status if condition is met.
#  define ARROW_RETURN_IF_(condition, status, expr)   \
    do {                                              \
      if (ARROW_PREDICT_FALSE(condition)) {           \
        ::arrow::Status _st = (status);               \
        _st.AddContextLine(__FILE__, __LINE__, expr); \
        return _st;                                   \
      }                                               \
    } while (0)

#else

#  define ARROW_RETURN_IF_(condition, status, _) \
    do {                                         \
      if (ARROW_PREDICT_FALSE(condition)) {      \
        return (status);                         \
      }                                          \
    } while (0)

#endif  // ARROW_EXTRA_ERROR_CONTEXT

#define ARROW_RETURN_IF(condition, status) \
  ARROW_RETURN_IF_(condition, status, ARROW_STRINGIFY(status))

/// \brief Propagate any non-successful Status to the caller
#define ARROW_RETURN_NOT_OK(status)                            \
  do {                                                         \
    ::arrow::Status __s = ::arrow::ToStatus(status);           \
    ARROW_RETURN_IF_(!__s.ok(), __s, ARROW_STRINGIFY(status)); \
  } while (false)

/// \brief Given `expr` and `warn_msg`; log `warn_msg` if `expr` is a non-ok status
#define ARROW_WARN_NOT_OK(expr, warn_msg) \
  do {                                    \
    ::arrow::Status _s = (expr);          \
    if (ARROW_PREDICT_FALSE(!_s.ok())) {  \
      _s.Warn(warn_msg);                  \
    }                                     \
  } while (false)

// This is an internal-use macro and should not be used in public headers.
#ifndef RETURN_NOT_OK
#  define RETURN_NOT_OK(s) ARROW_RETURN_NOT_OK(s)
#endif

namespace arrow {
namespace internal {

class StatusConstant;

}  // namespace internal

enum class StatusCode : char {
  OK = 0,
  OutOfMemory = 1,
  KeyError = 2,
  TypeError = 3,
  Invalid = 4,
  IOError = 5,
  CapacityError = 6,
  IndexError = 7,
  Cancelled = 8,
  UnknownError = 9,
  NotImplemented = 10,
  SerializationError = 11,
  RError = 13,
  // Gandiva range of errors
  CodeGenError = 40,
  ExpressionValidationError = 41,
  ExecutionError = 42,
  // Continue generic codes.
  AlreadyExists = 45
};

/// \brief An opaque class that allows subsystems to retain
/// additional information inside the Status.
class ARROW_EXPORT StatusDetail {
 public:
  virtual ~StatusDetail() = default;
  /// \brief Return a unique id for the type of the StatusDetail
  /// (effectively a poor man's substitute for RTTI).
  virtual const char* type_id() const = 0;
  /// \brief Produce a human-readable description of this status.
  virtual std::string ToString() const = 0;

  bool operator==(const StatusDetail& other) const noexcept {
    return std::string(type_id()) == other.type_id() && ToString() == other.ToString();
  }
};

/// \brief A type trait to declare a given type as Status-compatible.
///
/// This trait structure can be implemented if a type (such as Result<T>) embeds
/// error information that can be converted to the Status class.
/// It will make the given type usable directly in functions such as
/// Status::OrElse and error-checking macros such as ARROW_RETURN_NOT_OK.
template <typename T>
struct IntoStatus;

/// \brief Convert a Status-compatible object to Status
///
/// This generic function delegates to the IntoStatus type trait.
template <typename T>
constexpr decltype(auto) ToStatus(T&& t) {
  return IntoStatus<std::decay_t<T>>::ToStatus(std::forward<T>(t));
}

/// \brief Status outcome object (success or error)
///
/// The Status object is an object holding the outcome of an operation.
/// The outcome is represented as a StatusCode, either success
/// (StatusCode::OK) or an error (any other of the StatusCode enumeration values).
///
/// Additionally, if an error occurred, a specific error message is generally
/// attached.
class ARROW_EXPORT [[nodiscard]] Status : public util::EqualityComparable<Status>,
                                          public util::ToStringOstreamable<Status> {
 public:
  // Create a success status.
  constexpr Status() noexcept : state_(NULLPTR) {}
  ~Status() noexcept {
    if (ARROW_PREDICT_FALSE(state_ != NULL)) {
      if (!state_->is_constant) {
        DeleteState();
      }
    }
  }

  Status(StatusCode code, const std::string& msg);
  /// \brief Pluggable constructor for use by sub-systems.  detail cannot be null.
  Status(StatusCode code, std::string msg, std::shared_ptr<StatusDetail> detail);

  // Copy the specified status.
  inline Status(const Status& s);
  inline Status& operator=(const Status& s);

  // Move the specified status.
  inline Status(Status&& s) noexcept;
  inline Status& operator=(Status&& s) noexcept;

  inline bool Equals(const Status& s) const;

  // AND the statuses.
  inline Status operator&(const Status& s) const noexcept;
  inline Status operator&(Status&& s) const noexcept;
  inline Status& operator&=(const Status& s) noexcept;
  inline Status& operator&=(Status&& s) noexcept;

  /// Return a success status
  static Status OK() { return Status(); }

  template <typename... Args>
  static Status FromArgs(StatusCode code, Args&&... args) {
    return Status(code, internal::JoinToString(std::forward<Args>(args)...));
  }

  template <typename... Args>
  static Status FromDetailAndArgs(StatusCode code, std::shared_ptr<StatusDetail> detail,
                                  Args&&... args) {
    return Status(code, internal::JoinToString(std::forward<Args>(args)...),
                  std::move(detail));
  }

  /// Return an error status for out-of-memory conditions
  template <typename... Args>
  static Status OutOfMemory(Args&&... args) {
    return Status::FromArgs(StatusCode::OutOfMemory, std::forward<Args>(args)...);
  }

  /// Return an error status for failed key lookups (e.g. column name in a table)
  template <typename... Args>
  static Status KeyError(Args&&... args) {
    return Status::FromArgs(StatusCode::KeyError, std::forward<Args>(args)...);
  }

  /// Return an error status for type errors (such as mismatching data types)
  template <typename... Args>
  static Status TypeError(Args&&... args) {
    return Status::FromArgs(StatusCode::TypeError, std::forward<Args>(args)...);
  }

  /// Return an error status for unknown errors
  template <typename... Args>
  static Status UnknownError(Args&&... args) {
    return Status::FromArgs(StatusCode::UnknownError, std::forward<Args>(args)...);
  }

  /// Return an error status when an operation or a combination of operation and
  /// data types is unimplemented
  template <typename... Args>
  static Status NotImplemented(Args&&... args) {
    return Status::FromArgs(StatusCode::NotImplemented, std::forward<Args>(args)...);
  }

  /// Return an error status for invalid data (for example a string that fails parsing)
  template <typename... Args>
  static Status Invalid(Args&&... args) {
    return Status::FromArgs(StatusCode::Invalid, std::forward<Args>(args)...);
  }

  /// Return an error status for cancelled operation
  template <typename... Args>
  static Status Cancelled(Args&&... args) {
    return Status::FromArgs(StatusCode::Cancelled, std::forward<Args>(args)...);
  }

  /// Return an error status when an index is out of bounds
  template <typename... Args>
  static Status IndexError(Args&&... args) {
    return Status::FromArgs(StatusCode::IndexError, std::forward<Args>(args)...);
  }

  /// Return an error status when a container's capacity would exceed its limits
  template <typename... Args>
  static Status CapacityError(Args&&... args) {
    return Status::FromArgs(StatusCode::CapacityError, std::forward<Args>(args)...);
  }

  /// Return an error status when some IO-related operation failed
  template <typename... Args>
  static Status IOError(Args&&... args) {
    return Status::FromArgs(StatusCode::IOError, std::forward<Args>(args)...);
  }

  /// Return an error status when some (de)serialization operation failed
  template <typename... Args>
  static Status SerializationError(Args&&... args) {
    return Status::FromArgs(StatusCode::SerializationError, std::forward<Args>(args)...);
  }

  template <typename... Args>
  static Status RError(Args&&... args) {
    return Status::FromArgs(StatusCode::RError, std::forward<Args>(args)...);
  }

  template <typename... Args>
  static Status CodeGenError(Args&&... args) {
    return Status::FromArgs(StatusCode::CodeGenError, std::forward<Args>(args)...);
  }

  template <typename... Args>
  static Status ExpressionValidationError(Args&&... args) {
    return Status::FromArgs(StatusCode::ExpressionValidationError,
                            std::forward<Args>(args)...);
  }

  template <typename... Args>
  static Status ExecutionError(Args&&... args) {
    return Status::FromArgs(StatusCode::ExecutionError, std::forward<Args>(args)...);
  }

  template <typename... Args>
  static Status AlreadyExists(Args&&... args) {
    return Status::FromArgs(StatusCode::AlreadyExists, std::forward<Args>(args)...);
  }

  /// Return true iff the status indicates success.
  constexpr bool ok() const { return (state_ == NULLPTR); }

  /// Return true iff the status indicates an out-of-memory error.
  constexpr bool IsOutOfMemory() const { return code() == StatusCode::OutOfMemory; }
  /// Return true iff the status indicates a key lookup error.
  constexpr bool IsKeyError() const { return code() == StatusCode::KeyError; }
  /// Return true iff the status indicates invalid data.
  constexpr bool IsInvalid() const { return code() == StatusCode::Invalid; }
  /// Return true iff the status indicates a cancelled operation.
  constexpr bool IsCancelled() const { return code() == StatusCode::Cancelled; }
  /// Return true iff the status indicates an IO-related failure.
  constexpr bool IsIOError() const { return code() == StatusCode::IOError; }
  /// Return true iff the status indicates a container reaching capacity limits.
  constexpr bool IsCapacityError() const { return code() == StatusCode::CapacityError; }
  /// Return true iff the status indicates an out of bounds index.
  constexpr bool IsIndexError() const { return code() == StatusCode::IndexError; }
  /// Return true iff the status indicates a type error.
  constexpr bool IsTypeError() const { return code() == StatusCode::TypeError; }
  /// Return true iff the status indicates an unknown error.
  constexpr bool IsUnknownError() const { return code() == StatusCode::UnknownError; }
  /// Return true iff the status indicates an unimplemented operation.
  constexpr bool IsNotImplemented() const { return code() == StatusCode::NotImplemented; }
  /// Return true iff the status indicates a (de)serialization failure
  constexpr bool IsSerializationError() const {
    return code() == StatusCode::SerializationError;
  }
  /// Return true iff the status indicates a R-originated error.
  constexpr bool IsRError() const { return code() == StatusCode::RError; }

  constexpr bool IsCodeGenError() const { return code() == StatusCode::CodeGenError; }

  constexpr bool IsExpressionValidationError() const {
    return code() == StatusCode::ExpressionValidationError;
  }

  constexpr bool IsExecutionError() const { return code() == StatusCode::ExecutionError; }
  constexpr bool IsAlreadyExists() const { return code() == StatusCode::AlreadyExists; }

  /// \brief Return a string representation of this status suitable for printing.
  ///
  /// The string "OK" is returned for success.
  std::string ToString() const;

  /// \brief Return a string representation of this status without
  /// context lines suitable for printing.
  ///
  /// The string "OK" is returned for success.
  std::string ToStringWithoutContextLines() const;

  /// \brief Return a string representation of the status code, without the message
  /// text or POSIX code information.
  std::string CodeAsString() const;
  static std::string CodeAsString(StatusCode);

  /// \brief Return the StatusCode value attached to this status.
  constexpr StatusCode code() const { return ok() ? StatusCode::OK : state_->code; }

  /// \brief Return the specific error message attached to this status.
  const std::string& message() const;

  /// \brief Return the status detail attached to this message.
  const std::shared_ptr<StatusDetail>& detail() const;

  /// \brief Return a new Status copying the existing status, but
  /// updating with the existing detail.
  Status WithDetail(std::shared_ptr<StatusDetail> new_detail) const {
    return Status(code(), message(), std::move(new_detail));
  }

  /// \brief Return a new Status with changed message, copying the
  /// existing status code and detail.
  template <typename... Args>
  Status WithMessage(Args&&... args) const {
    return FromArgs(code(), std::forward<Args>(args)...).WithDetail(detail());
  }

  /// \brief Apply a functor if the status indicates an error
  ///
  /// This can be used to execute fallback or cleanup actions.
  ///
  /// If the status indicates a success, it is returned as-is.
  ///
  /// If the status indicates an error, the given functor is called with the status
  /// as argument.
  /// If the functor returns a new Status, it is returned.
  /// If the functor returns a Status-compatible object such as Result<T>, it is
  /// converted to Status and returned.
  /// If the functor returns void, the original Status is returned.
  template <typename OnError>
  Status OrElse(OnError&& on_error) {
    using RT = decltype(on_error(Status()));
    if (ARROW_PREDICT_TRUE(ok())) {
      return *this;
    }
    if constexpr (std::is_void_v<RT>) {
      on_error(*this);
      return *this;
    } else {
      return ToStatus(on_error(*this));
    }
  }

  void Warn() const;
  void Warn(const std::string& message) const;

  [[noreturn]] void Abort() const;
  [[noreturn]] void Abort(const std::string& message) const;

#ifdef ARROW_EXTRA_ERROR_CONTEXT
  void AddContextLine(const char* filename, int line, const char* expr);
#endif

 private:
  struct State {
    StatusCode code;
    bool is_constant;
    std::string msg;
    std::shared_ptr<StatusDetail> detail;
  };
  // OK status has a `NULL` state_.  Otherwise, `state_` points to
  // a `State` structure containing the error code and message(s)
  State* state_;

  void DeleteState() noexcept {
    // ARROW-2400: On certain compilers, splitting off the slow path improves
    // performance significantly.
    delete state_;
  }
  void CopyFrom(const Status& s);
  inline void MoveFrom(Status& s);

  friend class internal::StatusConstant;
};

void Status::MoveFrom(Status& s) {
  if (ARROW_PREDICT_FALSE(state_ != NULL)) {
    if (!state_->is_constant) {
      DeleteState();
    }
  }
  state_ = s.state_;
  s.state_ = NULLPTR;
}

Status::Status(const Status& s) : state_{NULLPTR} { CopyFrom(s); }

Status& Status::operator=(const Status& s) {
  // The following condition catches both aliasing (when this == &s),
  // and the common case where both s and *this are ok.
  if (state_ != s.state_) {
    CopyFrom(s);
  }
  return *this;
}

Status::Status(Status&& s) noexcept : state_(s.state_) { s.state_ = NULLPTR; }

Status& Status::operator=(Status&& s) noexcept {
  MoveFrom(s);
  return *this;
}

bool Status::Equals(const Status& s) const {
  if (state_ == s.state_) {
    return true;
  }

  if (ok() || s.ok()) {
    return false;
  }

  if (detail() != s.detail()) {
    if ((detail() && !s.detail()) || (!detail() && s.detail())) {
      return false;
    }
    return *detail() == *s.detail();
  }

  return code() == s.code() && message() == s.message();
}

/// \cond FALSE
// (note: emits warnings on Doxygen < 1.8.15,
//  see https://github.com/doxygen/doxygen/issues/6295)
Status Status::operator&(const Status& s) const noexcept {
  if (ok()) {
    return s;
  } else {
    return *this;
  }
}

Status Status::operator&(Status&& s) const noexcept {
  if (ok()) {
    return std::move(s);
  } else {
    return *this;
  }
}

Status& Status::operator&=(const Status& s) noexcept {
  if (ok() && !s.ok()) {
    CopyFrom(s);
  }
  return *this;
}

Status& Status::operator&=(Status&& s) noexcept {
  if (ok() && !s.ok()) {
    MoveFrom(s);
  }
  return *this;
}
/// \endcond

template <>
struct IntoStatus<Status> {
  static constexpr const Status& ToStatus(const Status& st) { return st; }
  static constexpr Status&& ToStatus(Status&& st) { return std::move(st); }
};

}  // namespace arrow
