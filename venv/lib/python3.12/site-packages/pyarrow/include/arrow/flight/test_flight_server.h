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

#include <memory>

#include "arrow/flight/server.h"
#include "arrow/flight/type_fwd.h"
#include "arrow/flight/visibility.h"
#include "arrow/status.h"

namespace arrow::flight {

class ARROW_FLIGHT_EXPORT TestFlightServer : public FlightServerBase {
 public:
  static std::unique_ptr<FlightServerBase> Make();

  Status ListFlights(const ServerCallContext& context, const Criteria* criteria,
                     std::unique_ptr<FlightListing>* listings) override;

  Status GetFlightInfo(const ServerCallContext& context, const FlightDescriptor& request,
                       std::unique_ptr<FlightInfo>* out) override;

  Status DoGet(const ServerCallContext& context, const Ticket& request,
               std::unique_ptr<FlightDataStream>* data_stream) override;

  Status DoPut(const ServerCallContext&, std::unique_ptr<FlightMessageReader> reader,
               std::unique_ptr<FlightMetadataWriter> writer) override;

  Status DoExchange(const ServerCallContext& context,
                    std::unique_ptr<FlightMessageReader> reader,
                    std::unique_ptr<FlightMessageWriter> writer) override;

  // A simple example - act like DoGet.
  Status RunExchangeGet(std::unique_ptr<FlightMessageReader> reader,
                        std::unique_ptr<FlightMessageWriter> writer);

  // A simple example - act like DoPut
  Status RunExchangePut(std::unique_ptr<FlightMessageReader> reader,
                        std::unique_ptr<FlightMessageWriter> writer);

  // Read some number of record batches from the client, send a
  // metadata message back with the count, then echo the batches back.
  Status RunExchangeCounter(std::unique_ptr<FlightMessageReader> reader,
                            std::unique_ptr<FlightMessageWriter> writer);

  // Read int64 batches from the client, each time sending back a
  // batch with a running sum of columns.
  Status RunExchangeTotal(std::unique_ptr<FlightMessageReader> reader,
                          std::unique_ptr<FlightMessageWriter> writer);

  // Echo the client's messages back.
  Status RunExchangeEcho(std::unique_ptr<FlightMessageReader> reader,
                         std::unique_ptr<FlightMessageWriter> writer);

  // Regression test for ARROW-13253
  Status RunExchangeLargeBatch(std::unique_ptr<FlightMessageReader>,
                               std::unique_ptr<FlightMessageWriter> writer);

  Status RunAction1(const Action& action, std::unique_ptr<ResultStream>* out);

  Status RunAction2(std::unique_ptr<ResultStream>* out);

  Status ListIncomingHeaders(const ServerCallContext& context, const Action& action,
                             std::unique_ptr<ResultStream>* out);

  Status DoAction(const ServerCallContext& context, const Action& action,
                  std::unique_ptr<ResultStream>* out) override;

  Status ListActions(const ServerCallContext& context,
                     std::vector<ActionType>* out) override;

  Status GetSchema(const ServerCallContext& context, const FlightDescriptor& request,
                   std::unique_ptr<SchemaResult>* schema) override;
};

}  // namespace arrow::flight
