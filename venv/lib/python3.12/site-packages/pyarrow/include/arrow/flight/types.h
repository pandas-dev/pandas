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

// Data structure for Flight RPC.

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

#include "arrow/buffer.h"
#include "arrow/flight/type_fwd.h"
#include "arrow/flight/visibility.h"
#include "arrow/ipc/options.h"
#include "arrow/ipc/writer.h"
#include "arrow/result.h"
#include "arrow/status.h"

namespace arrow {

class RecordBatch;
class Schema;
class Table;

namespace ipc {
class DictionaryMemo;
}  // namespace ipc

namespace util {
class Uri;
}  // namespace util

namespace flight {

ARROW_FLIGHT_EXPORT
extern const char* kSchemeGrpc;
ARROW_FLIGHT_EXPORT
extern const char* kSchemeGrpcTcp;
ARROW_FLIGHT_EXPORT
extern const char* kSchemeGrpcUnix;
ARROW_FLIGHT_EXPORT
extern const char* kSchemeGrpcTls;

class FlightClient;
class FlightServerBase;

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
using Timestamp =
    std::chrono::time_point<std::chrono::system_clock, std::chrono::nanoseconds>;

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
#  pragma warning(push)
#  pragma warning(disable : 4275)
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
#  pragma warning(pop)
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

namespace internal {

template <typename T>
struct remove_unique_ptr {
  using type = T;
};

template <typename T>
struct remove_unique_ptr<std::unique_ptr<T>> {
  using type = T;
};

// Base CRTP type
template <class T>
struct BaseType {
 protected:
  using SuperT = BaseType<T>;
  using SelfT = typename remove_unique_ptr<T>::type;

  const SelfT& self() const { return static_cast<const SelfT&>(*this); }
  SelfT& self() { return static_cast<SelfT&>(*this); }

 public:
  BaseType() = default;

  friend bool operator==(const SelfT& left, const SelfT& right) {
    return left.Equals(right);
  }
  friend bool operator!=(const SelfT& left, const SelfT& right) {
    return !left.Equals(right);
  }

  /// \brief Serialize this message to its wire-format representation.
  inline arrow::Result<std::string> SerializeToString() const {
    std::string out;
    ARROW_RETURN_NOT_OK(self().SelfT::SerializeToString(&out));
    return out;
  }

  inline static arrow::Result<T> Deserialize(std::string_view serialized) {
    T out;
    ARROW_RETURN_NOT_OK(SelfT::Deserialize(serialized, &out));
    return out;
  }

  inline arrow::Result<std::shared_ptr<Buffer>> SerializeToBuffer() const {
    std::string out;
    ARROW_RETURN_NOT_OK(self().SelfT::SerializeToString(&out));
    return Buffer::FromString(std::move(out));
  }
};

}  // namespace internal

//------------------------------------------------------------
// Wrapper types for Flight RPC protobuf messages

// A wrapper around arrow.flight.protocol.HandshakeRequest is not defined
// A wrapper around arrow.flight.protocol.HandshakeResponse is not defined

/// \brief message for simple auth
struct ARROW_FLIGHT_EXPORT BasicAuth : public internal::BaseType<BasicAuth> {
  std::string username;
  std::string password;

  BasicAuth() = default;
  BasicAuth(std::string username, std::string password)
      : username(std::move(username)), password(std::move(password)) {}

  std::string ToString() const;
  bool Equals(const BasicAuth& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized, BasicAuth* out);
};

// A wrapper around arrow.flight.protocol.Empty is not defined

/// \brief A type of action that can be performed with the DoAction RPC.
struct ARROW_FLIGHT_EXPORT ActionType : public internal::BaseType<ActionType> {
  /// \brief The name of the action.
  std::string type;

  /// \brief A human-readable description of the action.
  std::string description;

  ActionType() = default;

  ActionType(std::string type, std::string description)
      : type(std::move(type)), description(std::move(description)) {}

  std::string ToString() const;
  bool Equals(const ActionType& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized, ActionType* out);

  static const ActionType kCancelFlightInfo;
  static const ActionType kRenewFlightEndpoint;
  static const ActionType kSetSessionOptions;
  static const ActionType kGetSessionOptions;
  static const ActionType kCloseSession;
};

/// \brief Opaque selection criteria for ListFlights RPC
struct ARROW_FLIGHT_EXPORT Criteria : public internal::BaseType<Criteria> {
  /// Opaque criteria expression, dependent on server implementation
  std::string expression;

  Criteria() = default;
  Criteria(std::string expression)  // NOLINT runtime/explicit
      : expression(std::move(expression)) {}

  std::string ToString() const;
  bool Equals(const Criteria& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized, Criteria* out);
};

/// \brief An action to perform with the DoAction RPC
struct ARROW_FLIGHT_EXPORT Action : public internal::BaseType<Action> {
  /// The action type
  std::string type;

  /// The action content as a Buffer
  std::shared_ptr<Buffer> body;

  Action() = default;
  Action(std::string type, std::shared_ptr<Buffer> body)
      : type(std::move(type)), body(std::move(body)) {}

  std::string ToString() const;
  bool Equals(const Action& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized, Action* out);
};

/// \brief Opaque result returned after executing an action
struct ARROW_FLIGHT_EXPORT Result : public internal::BaseType<Result> {
  std::shared_ptr<Buffer> body;

  Result() = default;
  Result(std::shared_ptr<Buffer> body)  // NOLINT runtime/explicit
      : body(std::move(body)) {}

  std::string ToString() const;
  bool Equals(const Result& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized, Result* out);
};

/// \brief Schema result returned after a schema request RPC
struct ARROW_FLIGHT_EXPORT SchemaResult : public internal::BaseType<SchemaResult> {
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

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized, SchemaResult* out);

 private:
  std::string raw_schema_;
};

/// \brief A request to retrieve or generate a dataset
struct ARROW_FLIGHT_EXPORT FlightDescriptor
    : public internal::BaseType<FlightDescriptor> {
  enum DescriptorType {
    UNKNOWN = 0,  /// Unused
    PATH = 1,     /// Named path identifying a dataset
    CMD = 2       /// Opaque command to generate a dataset
  };

  /// The descriptor type
  DescriptorType type = UNKNOWN;

  /// Opaque value used to express a command. Should only be defined when type
  /// is CMD
  std::string cmd;

  /// List of strings identifying a particular dataset. Should only be defined
  /// when type is PATH
  std::vector<std::string> path;

  FlightDescriptor();
  FlightDescriptor(DescriptorType type, std::string cmd,
                   std::vector<std::string> path) noexcept;
  ~FlightDescriptor();

  /// \brief Get a human-readable form of this descriptor.
  std::string ToString() const;
  bool Equals(const FlightDescriptor& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Get the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Parse the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized, FlightDescriptor* out);

  // Convenience factory functions

  static FlightDescriptor Command(std::string cmd) {
    return FlightDescriptor{CMD, std::move(cmd), {}};
  }

  static FlightDescriptor Path(std::vector<std::string> path) {
    return FlightDescriptor{PATH, "", std::move(path)};
  }
};

/// \brief Data structure providing an opaque identifier or credential to use
/// when requesting a data stream with the DoGet RPC
struct ARROW_FLIGHT_EXPORT Ticket : public internal::BaseType<Ticket> {
  std::string ticket;

  Ticket() = default;
  Ticket(std::string ticket)  // NOLINT runtime/explicit
      : ticket(std::move(ticket)) {}

  std::string ToString() const;
  bool Equals(const Ticket& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Get the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Parse the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized, Ticket* out);
};

/// \brief A host location (a URI)
struct ARROW_FLIGHT_EXPORT Location : public internal::BaseType<Location> {
 public:
  /// \brief Initialize a blank location.
  Location();

  ~Location();

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

  /// \brief Get the scheme of this URI.
  std::string scheme() const;

  /// \brief Get a representation of this URI as a string.
  std::string ToString() const;
  bool Equals(const Location& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized, Location* out);

 private:
  friend class FlightClient;
  friend class FlightServerBase;
  std::shared_ptr<arrow::util::Uri> uri_;
};

/// \brief A flight ticket and list of locations where the ticket can be
/// redeemed
struct ARROW_FLIGHT_EXPORT FlightEndpoint : public internal::BaseType<FlightEndpoint> {
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

  FlightEndpoint() = default;
  FlightEndpoint(Ticket ticket, std::vector<Location> locations,
                 std::optional<Timestamp> expiration_time, std::string app_metadata)
      : ticket(std::move(ticket)),
        locations(std::move(locations)),
        expiration_time(expiration_time),
        app_metadata(std::move(app_metadata)) {}

  std::string ToString() const;
  bool Equals(const FlightEndpoint& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized, FlightEndpoint* out);
};

/// \brief The access coordinates for retrieval of a dataset, returned by
/// GetFlightInfo
class ARROW_FLIGHT_EXPORT FlightInfo
    : public internal::BaseType<std::unique_ptr<FlightInfo>> {
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

  /// \brief Factory method to construct a FlightInfo.
  static arrow::Result<FlightInfo> Make(const std::shared_ptr<Schema>& schema,
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
  /// \return Arrow result with the reconstructed Schema. Note that the schema
  ///   may be nullptr, as the schema is optional.
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

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Get the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Parse the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized,
                                   std::unique_ptr<FlightInfo>* out);

  std::string ToString() const;

  /// Compare two FlightInfo for equality. This will compare the
  /// serialized schema representations, NOT the logical equality of
  /// the schemas.
  bool Equals(const FlightInfo& other) const;

 private:
  Data data_;
  mutable std::shared_ptr<Schema> schema_;
  mutable bool reconstructed_schema_;
};

/// \brief The information to process a long-running query.
class ARROW_FLIGHT_EXPORT PollInfo
    : public internal::BaseType<std::unique_ptr<PollInfo>> {
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

  PollInfo(std::unique_ptr<FlightInfo> info, std::optional<FlightDescriptor> descriptor,
           std::optional<double> progress, std::optional<Timestamp> expiration_time)
      : info(std::move(info)),
        descriptor(std::move(descriptor)),
        progress(progress),
        expiration_time(expiration_time) {}

  PollInfo(const PollInfo& other)
      : info(other.info ? std::make_unique<FlightInfo>(*other.info) : NULLPTR),
        descriptor(other.descriptor),
        progress(other.progress),
        expiration_time(other.expiration_time) {}
  PollInfo(PollInfo&& other) noexcept = default;
  ~PollInfo() = default;
  PollInfo& operator=(const PollInfo& other) {
    info = other.info ? std::make_unique<FlightInfo>(*other.info) : NULLPTR;
    descriptor = other.descriptor;
    progress = other.progress;
    expiration_time = other.expiration_time;
    return *this;
  }
  PollInfo& operator=(PollInfo&& other) = default;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Get the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Parse the wire-format representation of this type.
  ///
  /// Useful when interoperating with non-Flight systems (e.g. REST
  /// services) that may want to return Flight types.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized,
                                   std::unique_ptr<PollInfo>* out);

  std::string ToString() const;

  /// Compare two PollInfo for equality. This will compare the
  /// serialized schema representations, NOT the logical equality of
  /// the schemas.
  bool Equals(const PollInfo& other) const;
};

/// \brief The request of the CancelFlightInfoRequest action.
struct ARROW_FLIGHT_EXPORT CancelFlightInfoRequest
    : public internal::BaseType<CancelFlightInfoRequest> {
  std::unique_ptr<FlightInfo> info;

  CancelFlightInfoRequest() = default;
  CancelFlightInfoRequest(std::unique_ptr<FlightInfo> info)  // NOLINT runtime/explicit
      : info(std::move(info)) {}

  std::string ToString() const;
  bool Equals(const CancelFlightInfoRequest& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized,
                                   CancelFlightInfoRequest* out);
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
struct ARROW_FLIGHT_EXPORT CancelFlightInfoResult
    : public internal::BaseType<CancelFlightInfoResult> {
  CancelStatus status = CancelStatus::kUnspecified;

  CancelFlightInfoResult() = default;
  CancelFlightInfoResult(CancelStatus status)  // NOLINT runtime/explicit
      : status(status) {}

  std::string ToString() const;
  bool Equals(const CancelFlightInfoResult& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized,
                                   CancelFlightInfoResult* out);
};

ARROW_FLIGHT_EXPORT
std::ostream& operator<<(std::ostream& os, CancelStatus status);

/// \brief The request of the RenewFlightEndpoint action.
struct ARROW_FLIGHT_EXPORT RenewFlightEndpointRequest
    : public internal::BaseType<RenewFlightEndpointRequest> {
  FlightEndpoint endpoint;

  RenewFlightEndpointRequest() = default;
  explicit RenewFlightEndpointRequest(FlightEndpoint endpoint)
      : endpoint(std::move(endpoint)) {}

  std::string ToString() const;
  bool Equals(const RenewFlightEndpointRequest& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized,
                                   RenewFlightEndpointRequest* out);
};

// FlightData in Flight.proto maps to FlightPayload here.

/// \brief Staging data structure for messages about to be put on the wire
///
/// This structure corresponds to FlightData in the protocol.
struct ARROW_FLIGHT_EXPORT FlightPayload {
  std::shared_ptr<Buffer> descriptor;
  std::shared_ptr<Buffer> app_metadata;
  ipc::IpcPayload ipc_message;

  FlightPayload() = default;
  FlightPayload(std::shared_ptr<Buffer> descriptor, std::shared_ptr<Buffer> app_metadata,
                ipc::IpcPayload ipc_message)
      : descriptor(std::move(descriptor)),
        app_metadata(std::move(app_metadata)),
        ipc_message(std::move(ipc_message)) {}

  /// \brief Check that the payload can be written to the wire.
  Status Validate() const;
};

// A wrapper around arrow.flight.protocol.PutResult is not defined

// Session management messages

/// \brief Variant supporting all possible value types for {Set,Get}SessionOptions
///
/// By convention, an attempt to set a valueless (std::monostate) SessionOptionValue
/// should attempt to unset or clear the named option value on the server.
using SessionOptionValue = std::variant<std::monostate, std::string, bool, int64_t,
                                        double, std::vector<std::string>>;
std::ostream& operator<<(std::ostream& os, const SessionOptionValue& v);

/// \brief A request to set a set of session options by name/value.
struct ARROW_FLIGHT_EXPORT SetSessionOptionsRequest
    : public internal::BaseType<SetSessionOptionsRequest> {
  std::map<std::string, SessionOptionValue> session_options;

  SetSessionOptionsRequest() = default;
  explicit SetSessionOptionsRequest(
      std::map<std::string, SessionOptionValue> session_options)
      : session_options(std::move(session_options)) {}

  std::string ToString() const;
  bool Equals(const SetSessionOptionsRequest& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized,
                                   SetSessionOptionsRequest* out);
};

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

/// \brief The result(s) of setting session option(s).
struct ARROW_FLIGHT_EXPORT SetSessionOptionsResult
    : public internal::BaseType<SetSessionOptionsResult> {
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

  SetSessionOptionsResult() = default;
  SetSessionOptionsResult(std::map<std::string, Error> errors)  // NOLINT runtime/explicit
      : errors(std::move(errors)) {}

  std::string ToString() const;
  bool Equals(const SetSessionOptionsResult& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized,
                                   SetSessionOptionsResult* out);
};

/// \brief A request to get current session options.
struct ARROW_FLIGHT_EXPORT GetSessionOptionsRequest
    : public internal::BaseType<GetSessionOptionsRequest> {
  GetSessionOptionsRequest() = default;

  std::string ToString() const;
  bool Equals(const GetSessionOptionsRequest& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized,
                                   GetSessionOptionsRequest* out);
};

/// \brief The current session options.
struct ARROW_FLIGHT_EXPORT GetSessionOptionsResult
    : public internal::BaseType<GetSessionOptionsResult> {
  std::map<std::string, SessionOptionValue> session_options;

  GetSessionOptionsResult() = default;
  GetSessionOptionsResult(  // NOLINT runtime/explicit
      std::map<std::string, SessionOptionValue> session_options)
      : session_options(std::move(session_options)) {}

  std::string ToString() const;
  bool Equals(const GetSessionOptionsResult& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized,
                                   GetSessionOptionsResult* out);
};

/// \brief A request to close the open client session.
struct ARROW_FLIGHT_EXPORT CloseSessionRequest
    : public internal::BaseType<CloseSessionRequest> {
  CloseSessionRequest() = default;

  std::string ToString() const;
  bool Equals(const CloseSessionRequest& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized, CloseSessionRequest* out);
};

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

/// \brief The result of attempting to close the client session.
struct ARROW_FLIGHT_EXPORT CloseSessionResult
    : public internal::BaseType<CloseSessionResult> {
  CloseSessionStatus status;

  CloseSessionResult() = default;
  CloseSessionResult(CloseSessionStatus status)  // NOLINT runtime/explicit
      : status(status) {}

  std::string ToString() const;
  bool Equals(const CloseSessionResult& other) const;

  using SuperT::Deserialize;
  using SuperT::SerializeToString;

  /// \brief Serialize this message to its wire-format representation.
  ///
  /// Use `SerializeToString()` if you want a Result-returning version.
  arrow::Status SerializeToString(std::string* out) const;

  /// \brief Deserialize this message from its wire-format representation.
  ///
  /// Use `Deserialize(serialized)` if you want a Result-returning version.
  static arrow::Status Deserialize(std::string_view serialized, CloseSessionResult* out);
};

//------------------------------------------------------------

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
  FlightStreamChunk() noexcept;
  ~FlightStreamChunk();

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
