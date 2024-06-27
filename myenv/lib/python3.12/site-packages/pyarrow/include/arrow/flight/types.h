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

// Data structure for Flight RPC. API should be considered experimental for now

#pragma once

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "arrow/flight/type_fwd.h"
#include "arrow/flight/visibility.h"
#include "arrow/ipc/options.h"
#include "arrow/ipc/writer.h"
#include "arrow/result.h"
#include "arrow/status.h"

namespace arrow {

class Buffer;
class RecordBatch;
class Schema;
class Status;
class Table;

namespace ipc {

class DictionaryMemo;

}  // namespace ipc

namespace util {

class Uri;

}  // namespace util

namespace flight {

/// \brief A timestamp compatible with Protocol Buffer's
/// google.protobuf.Timestamp:
///
/// https://protobuf.dev/reference/protobuf/google.protobuf/#timestamp
///
/// > A Timestamp represents a point in time independent of any time
/// > zone or calendar, represented as seconds and fractions of
/// > seconds at nanosecond resolution in UTC Epoch time. It is
/// > encoded using the Proleptic Gregorian Calendar which extends the
/// > Gregorian calendar backwards to year one. It is encoded assuming
/// > all minutes are 60 seconds long, i.e. leap seconds are "smeared"
/// > so that no leap second table is needed for interpretation. Range
/// > is from 0001-01-01T00:00:00Z to 9999-12-31T23:59:59.999999999Z.
using Timestamp = std::chrono::system_clock::time_point;

/// \brief A Flight-specific status code.  Used to encode some
///   additional status codes into an Arrow Status.
enum class FlightStatusCode : int8_t {
  /// An implementation error has occurred.
  Internal,
  /// A request timed out.
  TimedOut,
  /// A request was cancelled.
  Cancelled,
  /// We are not authenticated to the remote service.
  Unauthenticated,
  /// We do not have permission to make this request.
  Unauthorized,
  /// The remote service cannot handle this request at the moment.
  Unavailable,
  /// A request failed for some other reason
  Failed
};

// Silence warning
// "non dll-interface class RecordBatchReader used as base for dll-interface class"
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4275)
#endif

/// \brief Flight-specific error information in a Status.
class ARROW_FLIGHT_EXPORT FlightStatusDetail : public arrow::StatusDetail {
 public:
  explicit FlightStatusDetail(FlightStatusCode code) : code_{code} {}
  explicit FlightStatusDetail(FlightStatusCode code, std::string extra_info)
      : code_{code}, extra_info_(std::move(extra_info)) {}
  const char* type_id() const override;
  std::string ToString() const override;

  /// \brief Get the Flight status code.
  FlightStatusCode code() const;
  /// \brief Get the extra error info
  std::string extra_info() const;
  /// \brief Get the human-readable name of the status code.
  std::string CodeAsString() const;
  /// \brief Set the extra error info
  void set_extra_info(std::string extra_info);

  /// \brief Try to extract a \a FlightStatusDetail from any Arrow
  /// status.
  ///
  /// \return a \a FlightStatusDetail if it could be unwrapped, \a
  /// nullptr otherwise
  static std::shared_ptr<FlightStatusDetail> UnwrapStatus(const arrow::Status& status);

 private:
  FlightStatusCode code_;
  std::string extra_info_;
};

#ifdef _MSC_VER
#pragma warning(pop)
#endif

/// \brief Make an appropriate Arrow status for the given
/// Flight-specific status.
///
/// \param code The Flight status code.
/// \param message The message for the error.
/// \param extra_info Optional extra binary info for the error (eg protobuf)
ARROW_FLIGHT_EXPORT
Status MakeFlightError(FlightStatusCode code, std::string message,
                       std::string extra_info = {});

/// \brief Headers sent from the client or server.
///
/// Header values are ordered.
using CallHeaders = std::multimap<std::string_view, std::string_view>;

/// \brief A TLS certificate plus key.
struct ARROW_FLIGHT_EXPORT CertKeyPair {
  /// \brief The certificate in PEM format.
  std::string pem_cert;

  /// \brief The key in PEM format.
  std::string pem_key;
};

/// \brief A type of action that can be performed with the DoAction RPC.
struct ARROW_FLIGHT_EXPORT ActionType {
  /// \brief The name of the action.
  std::string type;

  /// \brief A human-readable description of the action.
  std::string description;

  std::string ToString() const;
  bool Equals(const ActionType& other) const;

  friend bool operator==(const ActionType& left, const ActionType& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const ActionType& left, const ActionType& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<ActionType> Deserialize(std::string_view serialized);

  static const ActionType kCancelFlightInfo;
  static const ActionType kRenewFlightEndpoint;
  static const ActionType kSetSessionOptions;
  static const ActionType kGetSessionOptions;
  static const ActionType kCloseSession;
};

/// \brief Opaque selection criteria for ListFlights RPC
struct ARROW_FLIGHT_EXPORT Criteria {
  /// Opaque criteria expression, dependent on server implementation
  std::string expression;

  std::string ToString() const;
  bool Equals(const Criteria& other) const;

  friend bool operator==(const Criteria& left, const Criteria& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const Criteria& left, const Criteria& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<Criteria> Deserialize(std::string_view serialized);
};

/// \brief An action to perform with the DoAction RPC
struct ARROW_FLIGHT_EXPORT Action {
  /// The action type
  std::string type;

  /// The action content as a Buffer
  std::shared_ptr<Buffer> body;

  std::string ToString() const;
  bool Equals(const Action& other) const;

  friend bool operator==(const Action& left, const Action& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const Action& left, const Action& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<Action> Deserialize(std::string_view serialized);
};

/// \brief Opaque result returned after executing an action
struct ARROW_FLIGHT_EXPORT Result {
  std::shared_ptr<Buffer> body;

  std::string ToString() const;
  bool Equals(const Result& other) const;

  friend bool operator==(const Result& left, const Result& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const Result& left, const Result& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<Result> Deserialize(std::string_view serialized);
};

enum class CancelStatus {
  /// The cancellation status is unknown. Servers should avoid using
  /// this value (send a kNotCancellable if the requested FlightInfo
  /// is not known). Clients can retry the request.
  kUnspecified = 0,
  /// The cancellation request is complete. Subsequent requests with
  /// the same payload may return kCancelled or a kNotCancellable error.
  kCancelled = 1,
  /// The cancellation request is in progress. The client may retry
  /// the cancellation request.
  kCancelling = 2,
  // The FlightInfo is not cancellable. The client should not retry the
  // cancellation request.
  kNotCancellable = 3,
};

/// \brief The result of the CancelFlightInfo action.
struct ARROW_FLIGHT_EXPORT CancelFlightInfoResult {
  CancelStatus status;

  std::string ToString() const;
  bool Equals(const CancelFlightInfoResult& other) const;

  friend bool operator==(const CancelFlightInfoResult& left,
                         const CancelFlightInfoResult& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const CancelFlightInfoResult& left,
                         const CancelFlightInfoResult& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<CancelFlightInfoResult> Deserialize(std::string_view serialized);
};

ARROW_FLIGHT_EXPORT
std::ostream& operator<<(std::ostream& os, CancelStatus status);

/// \brief message for simple auth
struct ARROW_FLIGHT_EXPORT BasicAuth {
  std::string username;
  std::string password;

  std::string ToString() const;
  bool Equals(const BasicAuth& other) const;

  friend bool operator==(const BasicAuth& left, const BasicAuth& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const BasicAuth& left, const BasicAuth& right) {
    return !(left == right);
  }

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<BasicAuth> Deserialize(std::string_view serialized);
  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;
};

/// \brief A request to retrieve or generate a dataset
struct ARROW_FLIGHT_EXPORT FlightDescriptor {
  enum DescriptorType {
    UNKNOWN = 0,  /// Unused
    PATH = 1,     /// Named path identifying a dataset
    CMD = 2       /// Opaque command to generate a dataset
  };

  /// The descriptor type
  DescriptorType type;

  /// Opaque value used to express a command. Should only be defined when type
  /// is CMD
  std::string cmd;

  /// List of strings identifying a particular dataset. Should only be defined
  /// when type is PATH
  std::vector<std::string> path;

  bool Equals(const FlightDescriptor& other) const;

  /// \brief Get a human-readable form of this descriptor.
  std::string ToString() const;

  /// \brief Get the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Parse the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  static arrow::Result<FlightDescriptor> Deserialize(std::string_view serialized);

  // Convenience factory functions

  static FlightDescriptor Command(const std::string& c) {
    return FlightDescriptor{CMD, c, {}};
  }

  static FlightDescriptor Path(const std::vector<std::string>& p) {
    return FlightDescriptor{PATH, "", p};
  }

  friend bool operator==(const FlightDescriptor& left, const FlightDescriptor& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const FlightDescriptor& left, const FlightDescriptor& right) {
    return !(left == right);
  }
};

/// \brief Data structure providing an opaque identifier or credential to use
/// when requesting a data stream with the DoGet RPC
struct ARROW_FLIGHT_EXPORT Ticket {
  std::string ticket;

  std::string ToString() const;
  bool Equals(const Ticket& other) const;

  friend bool operator==(const Ticket& left, const Ticket& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const Ticket& left, const Ticket& right) {
    return !(left == right);
  }

  /// \brief Get the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Parse the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  static arrow::Result<Ticket> Deserialize(std::string_view serialized);
};

class FlightClient;
class FlightServerBase;

ARROW_FLIGHT_EXPORT
extern const char* kSchemeGrpc;
ARROW_FLIGHT_EXPORT
extern const char* kSchemeGrpcTcp;
ARROW_FLIGHT_EXPORT
extern const char* kSchemeGrpcUnix;
ARROW_FLIGHT_EXPORT
extern const char* kSchemeGrpcTls;

/// \brief A host location (a URI)
struct ARROW_FLIGHT_EXPORT Location {
 public:
  /// \brief Initialize a blank location.
  Location();

  /// \brief Initialize a location by parsing a URI string
  static arrow::Result<Location> Parse(const std::string& uri_string);

  /// \brief Get the fallback URI.
  ///
  /// arrow-flight-reuse-connection://? means that a client may attempt to
  /// reuse an existing connection to a Flight service to fetch data instead
  /// of creating a new connection to one of the other locations listed in a
  /// FlightEndpoint response.
  static const Location& ReuseConnection();

  /// \brief Initialize a location for a non-TLS, gRPC-based Flight
  /// service from a host and port
  /// \param[in] host The hostname to connect to
  /// \param[in] port The port
  /// \return Arrow result with the resulting location
  static arrow::Result<Location> ForGrpcTcp(const std::string& host, const int port);

  /// \brief Initialize a location for a TLS-enabled, gRPC-based Flight
  /// service from a host and port
  /// \param[in] host The hostname to connect to
  /// \param[in] port The port
  /// \return Arrow result with the resulting location
  static arrow::Result<Location> ForGrpcTls(const std::string& host, const int port);

  /// \brief Initialize a location for a domain socket-based Flight
  /// service
  /// \param[in] path The path to the domain socket
  /// \return Arrow result with the resulting location
  static arrow::Result<Location> ForGrpcUnix(const std::string& path);

  /// \brief Initialize a location based on a URI scheme
  static arrow::Result<Location> ForScheme(const std::string& scheme,
                                           const std::string& host, const int port);

  /// \brief Get a representation of this URI as a string.
  std::string ToString() const;

  /// \brief Get the scheme of this URI.
  std::string scheme() const;

  bool Equals(const Location& other) const;

  friend bool operator==(const Location& left, const Location& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const Location& left, const Location& right) {
    return !(left == right);
  }

 private:
  friend class FlightClient;
  friend class FlightServerBase;
  std::shared_ptr<arrow::util::Uri> uri_;
};

/// \brief A flight ticket and list of locations where the ticket can be
/// redeemed
struct ARROW_FLIGHT_EXPORT FlightEndpoint {
  /// Opaque ticket identify; use with DoGet RPC
  Ticket ticket;

  /// List of locations where ticket can be redeemed. If the list is empty, the
  /// ticket can only be redeemed on the current service where the ticket was
  /// generated
  std::vector<Location> locations;

  /// Expiration time of this stream. If present, clients may assume
  /// they can retry DoGet requests. Otherwise, clients should avoid
  /// retrying DoGet requests.
  std::optional<Timestamp> expiration_time;

  /// Opaque Application-defined metadata
  std::string app_metadata;

  std::string ToString() const;
  bool Equals(const FlightEndpoint& other) const;

  friend bool operator==(const FlightEndpoint& left, const FlightEndpoint& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const FlightEndpoint& left, const FlightEndpoint& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<FlightEndpoint> Deserialize(std::string_view serialized);
};

/// \brief The request of the RenewFlightEndpoint action.
struct ARROW_FLIGHT_EXPORT RenewFlightEndpointRequest {
  FlightEndpoint endpoint;

  std::string ToString() const;
  bool Equals(const RenewFlightEndpointRequest& other) const;

  friend bool operator==(const RenewFlightEndpointRequest& left,
                         const RenewFlightEndpointRequest& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const RenewFlightEndpointRequest& left,
                         const RenewFlightEndpointRequest& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<RenewFlightEndpointRequest> Deserialize(
      std::string_view serialized);
};

/// \brief Staging data structure for messages about to be put on the wire
///
/// This structure corresponds to FlightData in the protocol.
struct ARROW_FLIGHT_EXPORT FlightPayload {
  std::shared_ptr<Buffer> descriptor;
  std::shared_ptr<Buffer> app_metadata;
  ipc::IpcPayload ipc_message;

  /// \brief Check that the payload can be written to the wire.
  Status Validate() const;
};

/// \brief Schema result returned after a schema request RPC
struct ARROW_FLIGHT_EXPORT SchemaResult {
 public:
  SchemaResult() = default;
  explicit SchemaResult(std::string schema) : raw_schema_(std::move(schema)) {}

  /// \brief Factory method to construct a SchemaResult.
  static arrow::Result<std::unique_ptr<SchemaResult>> Make(const Schema& schema);

  /// \brief return schema
  /// \param[in,out] dictionary_memo for dictionary bookkeeping, will
  /// be modified
  /// \return Arrow result with the reconstructed Schema
  arrow::Result<std::shared_ptr<Schema>> GetSchema(
      ipc::DictionaryMemo* dictionary_memo) const;

  const std::string& serialized_schema() const { return raw_schema_; }

  std::string ToString() const;
  bool Equals(const SchemaResult& other) const;

  friend bool operator==(const SchemaResult& left, const SchemaResult& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const SchemaResult& left, const SchemaResult& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<SchemaResult> Deserialize(std::string_view serialized);

 private:
  std::string raw_schema_;
};

/// \brief The access coordinates for retrieval of a dataset, returned by
/// GetFlightInfo
class ARROW_FLIGHT_EXPORT FlightInfo {
 public:
  struct Data {
    std::string schema;
    FlightDescriptor descriptor;
    std::vector<FlightEndpoint> endpoints;
    int64_t total_records = -1;
    int64_t total_bytes = -1;
    bool ordered = false;
    std::string app_metadata;
  };

  explicit FlightInfo(Data data) : data_(std::move(data)), reconstructed_schema_(false) {}

  /// \brief Factory method to construct a FlightInfo.
  static arrow::Result<FlightInfo> Make(const Schema& schema,
                                        const FlightDescriptor& descriptor,
                                        const std::vector<FlightEndpoint>& endpoints,
                                        int64_t total_records, int64_t total_bytes,
                                        bool ordered = false,
                                        std::string app_metadata = "");

  /// \brief Deserialize the Arrow schema of the dataset. Populate any
  ///   dictionary encoded fields into a DictionaryMemo for
  ///   bookkeeping
  /// \param[in,out] dictionary_memo for dictionary bookkeeping, will
  /// be modified
  /// \return Arrow result with the reconstructed Schema
  arrow::Result<std::shared_ptr<Schema>> GetSchema(
      ipc::DictionaryMemo* dictionary_memo) const;

  const std::string& serialized_schema() const { return data_.schema; }

  /// The descriptor associated with this flight, may not be set
  const FlightDescriptor& descriptor() const { return data_.descriptor; }

  /// A list of endpoints associated with the flight (dataset). To consume the
  /// whole flight, all endpoints must be consumed
  const std::vector<FlightEndpoint>& endpoints() const { return data_.endpoints; }

  /// The total number of records (rows) in the dataset. If unknown, set to -1
  int64_t total_records() const { return data_.total_records; }

  /// The total number of bytes in the dataset. If unknown, set to -1
  int64_t total_bytes() const { return data_.total_bytes; }

  /// Whether endpoints are in the same order as the data.
  bool ordered() const { return data_.ordered; }

  /// Application-defined opaque metadata
  const std::string& app_metadata() const { return data_.app_metadata; }

  /// \brief Get the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Parse the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  static arrow::Result<std::unique_ptr<FlightInfo>> Deserialize(
      std::string_view serialized);

  std::string ToString() const;

  /// Compare two FlightInfo for equality. This will compare the
  /// serialized schema representations, NOT the logical equality of
  /// the schemas.
  bool Equals(const FlightInfo& other) const;

  friend bool operator==(const FlightInfo& left, const FlightInfo& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const FlightInfo& left, const FlightInfo& right) {
    return !(left == right);
  }

 private:
  Data data_;
  mutable std::shared_ptr<Schema> schema_;
  mutable bool reconstructed_schema_;
};

/// \brief The information to process a long-running query.
class ARROW_FLIGHT_EXPORT PollInfo {
 public:
  /// The currently available results so far.
  std::unique_ptr<FlightInfo> info = NULLPTR;
  /// The descriptor the client should use on the next try. If unset,
  /// the query is complete.
  std::optional<FlightDescriptor> descriptor = std::nullopt;
  /// Query progress. Must be in [0.0, 1.0] but need not be
  /// monotonic or nondecreasing. If unknown, do not set.
  std::optional<double> progress = std::nullopt;
  /// Expiration time for this request. After this passes, the server
  /// might not accept the poll descriptor anymore (and the query may
  /// be cancelled). This may be updated on a call to PollFlightInfo.
  std::optional<Timestamp> expiration_time = std::nullopt;

  PollInfo()
      : info(NULLPTR),
        descriptor(std::nullopt),
        progress(std::nullopt),
        expiration_time(std::nullopt) {}

  explicit PollInfo(std::unique_ptr<FlightInfo> info,
                    std::optional<FlightDescriptor> descriptor,
                    std::optional<double> progress,
                    std::optional<Timestamp> expiration_time)
      : info(std::move(info)),
        descriptor(std::move(descriptor)),
        progress(progress),
        expiration_time(expiration_time) {}

  // Must not be explicit; to declare one we must declare all ("rule of five")
  PollInfo(const PollInfo& other)  // NOLINT(runtime/explicit)
      : info(other.info ? std::make_unique<FlightInfo>(*other.info) : NULLPTR),
        descriptor(other.descriptor),
        progress(other.progress),
        expiration_time(other.expiration_time) {}
  PollInfo(PollInfo&& other) noexcept = default;  // NOLINT(runtime/explicit)
  ~PollInfo() = default;
  PollInfo& operator=(const PollInfo& other) {
    info = other.info ? std::make_unique<FlightInfo>(*other.info) : NULLPTR;
    descriptor = other.descriptor;
    progress = other.progress;
    expiration_time = other.expiration_time;
    return *this;
  }
  PollInfo& operator=(PollInfo&& other) = default;

  /// \brief Get the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Parse the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  static arrow::Result<std::unique_ptr<PollInfo>> Deserialize(
      std::string_view serialized);

  std::string ToString() const;

  /// Compare two PollInfo for equality. This will compare the
  /// serialized schema representations, NOT the logical equality of
  /// the schemas.
  bool Equals(const PollInfo& other) const;

  friend bool operator==(const PollInfo& left, const PollInfo& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const PollInfo& left, const PollInfo& right) {
    return !(left == right);
  }
};

/// \brief The request of the CancelFlightInfoRequest action.
struct ARROW_FLIGHT_EXPORT CancelFlightInfoRequest {
  std::unique_ptr<FlightInfo> info;

  std::string ToString() const;
  bool Equals(const CancelFlightInfoRequest& other) const;

  friend bool operator==(const CancelFlightInfoRequest& left,
                         const CancelFlightInfoRequest& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const CancelFlightInfoRequest& left,
                         const CancelFlightInfoRequest& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<CancelFlightInfoRequest> Deserialize(std::string_view serialized);
};

/// \brief Variant supporting all possible value types for {Set,Get}SessionOptions
///
/// By convention, an attempt to set a valueless (std::monostate) SessionOptionValue
/// should attempt to unset or clear the named option value on the server.
using SessionOptionValue = std::variant<std::monostate, std::string, bool, int64_t,
                                        double, std::vector<std::string>>;

/// \brief The result of setting a session option.
enum class SetSessionOptionErrorValue : int8_t {
  /// \brief The status of setting the option is unknown.
  ///
  /// Servers should avoid using this value (send a NOT_FOUND error if the requested
  /// session is not known). Clients can retry the request.
  kUnspecified,
  /// \brief The given session option name is invalid.
  kInvalidName,
  /// \brief The session option value or type is invalid.
  kInvalidValue,
  /// \brief The session option cannot be set.
  kError
};
std::string ToString(const SetSessionOptionErrorValue& error_value);
std::ostream& operator<<(std::ostream& os, const SetSessionOptionErrorValue& error_value);

/// \brief The result of closing a session.
enum class CloseSessionStatus : int8_t {
  // \brief The session close status is unknown.
  //
  // Servers should avoid using this value (send a NOT_FOUND error if the requested
  // session is not known). Clients can retry the request.
  kUnspecified,
  // \brief The session close request is complete.
  //
  // Subsequent requests with the same session produce a NOT_FOUND error.
  kClosed,
  // \brief The session close request is in progress.
  //
  // The client may retry the request.
  kClosing,
  // \brief The session is not closeable.
  //
  // The client should not retry the request.
  kNotClosable
};
std::string ToString(const CloseSessionStatus& status);
std::ostream& operator<<(std::ostream& os, const CloseSessionStatus& status);

/// \brief A request to set a set of session options by name/value.
struct ARROW_FLIGHT_EXPORT SetSessionOptionsRequest {
  std::map<std::string, SessionOptionValue> session_options;

  std::string ToString() const;
  bool Equals(const SetSessionOptionsRequest& other) const;

  friend bool operator==(const SetSessionOptionsRequest& left,
                         const SetSessionOptionsRequest& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const SetSessionOptionsRequest& left,
                         const SetSessionOptionsRequest& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<SetSessionOptionsRequest> Deserialize(std::string_view serialized);
};

/// \brief The result(s) of setting session option(s).
struct ARROW_FLIGHT_EXPORT SetSessionOptionsResult {
  struct Error {
    SetSessionOptionErrorValue value;

    bool Equals(const Error& other) const { return value == other.value; }
    friend bool operator==(const Error& left, const Error& right) {
      return left.Equals(right);
    }
    friend bool operator!=(const Error& left, const Error& right) {
      return !(left == right);
    }
  };

  std::map<std::string, Error> errors;

  std::string ToString() const;
  bool Equals(const SetSessionOptionsResult& other) const;

  friend bool operator==(const SetSessionOptionsResult& left,
                         const SetSessionOptionsResult& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const SetSessionOptionsResult& left,
                         const SetSessionOptionsResult& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<SetSessionOptionsResult> Deserialize(std::string_view serialized);
};

/// \brief A request to get current session options.
struct ARROW_FLIGHT_EXPORT GetSessionOptionsRequest {
  std::string ToString() const;
  bool Equals(const GetSessionOptionsRequest& other) const;

  friend bool operator==(const GetSessionOptionsRequest& left,
                         const GetSessionOptionsRequest& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const GetSessionOptionsRequest& left,
                         const GetSessionOptionsRequest& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<GetSessionOptionsRequest> Deserialize(std::string_view serialized);
};

/// \brief The current session options.
struct ARROW_FLIGHT_EXPORT GetSessionOptionsResult {
  std::map<std::string, SessionOptionValue> session_options;

  std::string ToString() const;
  bool Equals(const GetSessionOptionsResult& other) const;

  friend bool operator==(const GetSessionOptionsResult& left,
                         const GetSessionOptionsResult& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const GetSessionOptionsResult& left,
                         const GetSessionOptionsResult& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<GetSessionOptionsResult> Deserialize(std::string_view serialized);
};

/// \brief A request to close the open client session.
struct ARROW_FLIGHT_EXPORT CloseSessionRequest {
  std::string ToString() const;
  bool Equals(const CloseSessionRequest& other) const;

  friend bool operator==(const CloseSessionRequest& left,
                         const CloseSessionRequest& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const CloseSessionRequest& left,
                         const CloseSessionRequest& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<CloseSessionRequest> Deserialize(std::string_view serialized);
};

/// \brief The result of attempting to close the client session.
struct ARROW_FLIGHT_EXPORT CloseSessionResult {
  CloseSessionStatus status;

  std::string ToString() const;
  bool Equals(const CloseSessionResult& other) const;

  friend bool operator==(const CloseSessionResult& left,
                         const CloseSessionResult& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const CloseSessionResult& left,
                         const CloseSessionResult& right) {
    return !(left == right);
  }

  /// \brief Serialize this message to its wire-format representation.
  arrow::Result<std::string> SerializeToString() const;

  /// \brief Deserialize this message from its wire-format representation.
  static arrow::Result<CloseSessionResult> Deserialize(std::string_view serialized);
};

/// \brief An iterator to FlightInfo instances returned by ListFlights.
class ARROW_FLIGHT_EXPORT FlightListing {
 public:
  virtual ~FlightListing() = default;

  /// \brief Retrieve the next FlightInfo from the iterator.
  /// \return Arrow result with a single FlightInfo. Set to \a nullptr if there
  /// are none left.
  virtual arrow::Result<std::unique_ptr<FlightInfo>> Next() = 0;
};

/// \brief An iterator to Result instances returned by DoAction.
class ARROW_FLIGHT_EXPORT ResultStream {
 public:
  virtual ~ResultStream() = default;

  /// \brief Retrieve the next Result from the iterator.
  /// \return Arrow result with a single Result. Set to \a nullptr if there are none left.
  virtual arrow::Result<std::unique_ptr<Result>> Next() = 0;

  /// \brief Read and drop the remaining messages to get the error (if any) from a server.
  /// \return Status OK if this is no error from a server, any other status if a
  /// server returns an error.
  Status Drain();
};

/// \brief A holder for a RecordBatch with associated Flight metadata.
struct ARROW_FLIGHT_EXPORT FlightStreamChunk {
 public:
  std::shared_ptr<RecordBatch> data;
  std::shared_ptr<Buffer> app_metadata;
};

/// \brief An interface to read Flight data with metadata.
class ARROW_FLIGHT_EXPORT MetadataRecordBatchReader {
 public:
  virtual ~MetadataRecordBatchReader() = default;

  /// \brief Get the schema for this stream.
  virtual arrow::Result<std::shared_ptr<Schema>> GetSchema() = 0;

  /// \brief Get the next message from Flight. If the stream is
  /// finished, then the members of \a FlightStreamChunk will be
  /// nullptr.
  virtual arrow::Result<FlightStreamChunk> Next() = 0;

  /// \brief Consume entire stream as a vector of record batches
  virtual arrow::Result<std::vector<std::shared_ptr<RecordBatch>>> ToRecordBatches();

  /// \brief Consume entire stream as a Table
  virtual arrow::Result<std::shared_ptr<Table>> ToTable();
};

/// \brief Convert a MetadataRecordBatchReader to a regular RecordBatchReader.
ARROW_FLIGHT_EXPORT
arrow::Result<std::shared_ptr<RecordBatchReader>> MakeRecordBatchReader(
    std::shared_ptr<MetadataRecordBatchReader> reader);

/// \brief An interface to write IPC payloads with metadata.
class ARROW_FLIGHT_EXPORT MetadataRecordBatchWriter : public ipc::RecordBatchWriter {
 public:
  virtual ~MetadataRecordBatchWriter() = default;
  /// \brief Begin writing data with the given schema. Only used with \a DoExchange.
  virtual Status Begin(const std::shared_ptr<Schema>& schema,
                       const ipc::IpcWriteOptions& options) = 0;
  virtual Status Begin(const std::shared_ptr<Schema>& schema);
  virtual Status WriteMetadata(std::shared_ptr<Buffer> app_metadata) = 0;
  virtual Status WriteWithMetadata(const RecordBatch& batch,
                                   std::shared_ptr<Buffer> app_metadata) = 0;
};

/// \brief A FlightListing implementation based on a vector of
/// FlightInfo objects.
///
/// This can be iterated once, then it is consumed.
class ARROW_FLIGHT_EXPORT SimpleFlightListing : public FlightListing {
 public:
  explicit SimpleFlightListing(const std::vector<FlightInfo>& flights);
  explicit SimpleFlightListing(std::vector<FlightInfo>&& flights);

  arrow::Result<std::unique_ptr<FlightInfo>> Next() override;

 private:
  int position_;
  std::vector<FlightInfo> flights_;
};

/// \brief A ResultStream implementation based on a vector of
/// Result objects.
///
/// This can be iterated once, then it is consumed.
class ARROW_FLIGHT_EXPORT SimpleResultStream : public ResultStream {
 public:
  explicit SimpleResultStream(std::vector<Result>&& results);
  arrow::Result<std::unique_ptr<Result>> Next() override;

 private:
  std::vector<Result> results_;
  size_t position_;
};

/// \defgroup flight-error Error Handling
/// Types for handling errors from RPCs.  Flight uses a set of status
/// codes standardized across Flight implementations, so these types
/// let applications work directly with those codes instead of having
/// to translate to and from Arrow Status.
/// @{

/// \brief Abstract status code for an RPC as per the Flight
///   specification.
enum class TransportStatusCode {
  /// \brief No error.
  kOk = 0,
  /// \brief An unknown error occurred.
  kUnknown = 1,
  /// \brief An error occurred in the transport implementation, or an
  ///   error internal to the service implementation occurred.
  kInternal = 2,
  /// \brief An argument is invalid.
  kInvalidArgument = 3,
  /// \brief The request timed out.
  kTimedOut = 4,
  /// \brief An argument is not necessarily invalid, but references
  ///   some resource that does not exist.  Prefer over
  ///   kInvalidArgument where applicable.
  kNotFound = 5,
  /// \brief The request attempted to create some resource that does
  ///   not exist.
  kAlreadyExists = 6,
  /// \brief The request was explicitly cancelled.
  kCancelled = 7,
  /// \brief The client is not authenticated.
  kUnauthenticated = 8,
  /// \brief The client is not authorized to perform this request.
  kUnauthorized = 9,
  /// \brief The request is not implemented
  kUnimplemented = 10,
  /// \brief There is a network connectivity error, or some resource
  ///   is otherwise unavailable.  Most likely a temporary condition.
  kUnavailable = 11,
};

/// \brief Convert a code to a string.
std::string ToString(TransportStatusCode code);

/// \brief An error from an RPC call, using Flight error codes directly
///   instead of trying to translate to Arrow Status.
///
/// Currently, only attached to the Status passed to AsyncListener::OnFinish.
///
/// This API is EXPERIMENTAL.
class ARROW_FLIGHT_EXPORT TransportStatusDetail : public StatusDetail {
 public:
  constexpr static const char* kTypeId = "flight::TransportStatusDetail";
  explicit TransportStatusDetail(TransportStatusCode code, std::string message,
                                 std::vector<std::pair<std::string, std::string>> details)
      : code_(code), message_(std::move(message)), details_(std::move(details)) {}
  const char* type_id() const override { return kTypeId; }
  std::string ToString() const override;

  static std::optional<std::reference_wrapper<const TransportStatusDetail>> Unwrap(
      const Status& status);

  TransportStatusCode code() const { return code_; }
  std::string_view message() const { return message_; }
  const std::vector<std::pair<std::string, std::string>>& details() const {
    return details_;
  }

 private:
  TransportStatusCode code_;
  std::string message_;
  std::vector<std::pair<std::string, std::string>> details_;
};

/// @}

}  // namespace flight
}  // namespace arrow
