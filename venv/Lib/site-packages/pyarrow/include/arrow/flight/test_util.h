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

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstdint>
#include <memory>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "arrow/status.h"
#include "arrow/testing/gtest_util.h"
#include "arrow/testing/process.h"
#include "arrow/testing/util.h"

#include "arrow/flight/client.h"
#include "arrow/flight/server.h"
#include "arrow/flight/types.h"
#include "arrow/flight/visibility.h"

namespace arrow {
namespace flight {

// ----------------------------------------------------------------------
// Helpers to compare values for equality

inline void AssertEqual(const FlightInfo& expected, const FlightInfo& actual) {
  ipc::DictionaryMemo expected_memo;
  ipc::DictionaryMemo actual_memo;
  ASSERT_OK_AND_ASSIGN(auto ex_schema, expected.GetSchema(&expected_memo));
  ASSERT_OK_AND_ASSIGN(auto actual_schema, actual.GetSchema(&actual_memo));

  AssertSchemaEqual(*ex_schema, *actual_schema);
  ASSERT_EQ(expected.total_records(), actual.total_records());
  ASSERT_EQ(expected.total_bytes(), actual.total_bytes());

  ASSERT_EQ(expected.descriptor(), actual.descriptor());
  ASSERT_THAT(actual.endpoints(), ::testing::ContainerEq(expected.endpoints()));
}

// ----------------------------------------------------------------------
// Fixture to use for running test servers

class ARROW_FLIGHT_EXPORT TestServer {
 public:
  explicit TestServer(const std::string& executable_name)
      : executable_name_(executable_name), port_(::arrow::GetListenPort()) {}
  TestServer(const std::string& executable_name, int port)
      : executable_name_(executable_name), port_(port) {}
  TestServer(const std::string& executable_name, const std::string& unix_sock)
      : executable_name_(executable_name), unix_sock_(unix_sock) {}

  Status Start(const std::vector<std::string>& extra_args);
  Status Start() { return Start({}); }

  void Stop();

  bool IsRunning();

  int port() const;
  const std::string& unix_sock() const;

 private:
  std::string executable_name_;
  int port_;
  std::string unix_sock_;
  std::unique_ptr<util::Process> server_process_;
};

// Helper to initialize a server and matching client with callbacks to
// populate options.
template <typename T, typename... Args>
Status MakeServer(const Location& location, std::unique_ptr<FlightServerBase>* server,
                  std::unique_ptr<FlightClient>* client,
                  std::function<Status(FlightServerOptions*)> make_server_options,
                  std::function<Status(FlightClientOptions*)> make_client_options,
                  Args&&... server_args) {
  *server = std::make_unique<T>(std::forward<Args>(server_args)...);
  FlightServerOptions server_options(location);
  RETURN_NOT_OK(make_server_options(&server_options));
  RETURN_NOT_OK((*server)->Init(server_options));
  std::string uri =
      location.scheme() + "://127.0.0.1:" + std::to_string((*server)->port());
  ARROW_ASSIGN_OR_RAISE(auto real_location, Location::Parse(uri));
  FlightClientOptions client_options = FlightClientOptions::Defaults();
  RETURN_NOT_OK(make_client_options(&client_options));
  return FlightClient::Connect(real_location, client_options).Value(client);
}

// Helper to initialize a server and matching client with callbacks to
// populate options.
template <typename T, typename... Args>
Status MakeServer(std::unique_ptr<FlightServerBase>* server,
                  std::unique_ptr<FlightClient>* client,
                  std::function<Status(FlightServerOptions*)> make_server_options,
                  std::function<Status(FlightClientOptions*)> make_client_options,
                  Args&&... server_args) {
  ARROW_ASSIGN_OR_RAISE(auto location, Location::ForGrpcTcp("localhost", 0));
  return MakeServer<T>(location, server, client, std::move(make_server_options),
                       std::move(make_client_options),
                       std::forward<Args>(server_args)...);
}

// ----------------------------------------------------------------------
// A FlightDataStream that numbers the record batches
/// \brief A basic implementation of FlightDataStream that will provide
/// a sequence of FlightData messages to be written to a stream
class ARROW_FLIGHT_EXPORT NumberingStream : public FlightDataStream {
 public:
  explicit NumberingStream(std::unique_ptr<FlightDataStream> stream);

  std::shared_ptr<Schema> schema() override;
  arrow::Result<FlightPayload> GetSchemaPayload() override;
  arrow::Result<FlightPayload> Next() override;

 private:
  int counter_;
  std::shared_ptr<FlightDataStream> stream_;
};

// ----------------------------------------------------------------------
// Example data for test-server and unit tests

ARROW_FLIGHT_EXPORT
std::shared_ptr<Schema> ExampleIntSchema();

ARROW_FLIGHT_EXPORT
std::shared_ptr<Schema> ExampleStringSchema();

ARROW_FLIGHT_EXPORT
std::shared_ptr<Schema> ExampleDictSchema();

ARROW_FLIGHT_EXPORT
std::shared_ptr<Schema> ExampleLargeSchema();

ARROW_FLIGHT_EXPORT
Status ExampleIntBatches(RecordBatchVector* out);

ARROW_FLIGHT_EXPORT
Status ExampleFloatBatches(RecordBatchVector* out);

ARROW_FLIGHT_EXPORT
Status ExampleDictBatches(RecordBatchVector* out);

ARROW_FLIGHT_EXPORT
Status ExampleNestedBatches(RecordBatchVector* out);

ARROW_FLIGHT_EXPORT
Status ExampleLargeBatches(RecordBatchVector* out);

ARROW_FLIGHT_EXPORT
arrow::Result<std::shared_ptr<RecordBatch>> VeryLargeBatch();

ARROW_FLIGHT_EXPORT
std::vector<FlightInfo> ExampleFlightInfo();

ARROW_FLIGHT_EXPORT
std::vector<ActionType> ExampleActionTypes();

ARROW_FLIGHT_EXPORT
FlightInfo MakeFlightInfo(const Schema& schema, const FlightDescriptor& descriptor,
                          const std::vector<FlightEndpoint>& endpoints,
                          int64_t total_records, int64_t total_bytes, bool ordered,
                          std::string app_metadata);

ARROW_FLIGHT_EXPORT
Status ExampleTlsCertificates(std::vector<CertKeyPair>* out);

ARROW_FLIGHT_EXPORT
Status ExampleTlsCertificateRoot(CertKeyPair* out);

}  // namespace flight
}  // namespace arrow
