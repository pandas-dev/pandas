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

#include <string>

#include "arrow/flight/client_auth.h"
#include "arrow/flight/server.h"
#include "arrow/flight/server_auth.h"
#include "arrow/flight/types.h"
#include "arrow/flight/visibility.h"
#include "arrow/status.h"

// A pair of authentication handlers that check for a predefined password
// and set the peer identity to a predefined username.

namespace arrow::flight {

class ARROW_FLIGHT_EXPORT TestServerAuthHandler : public ServerAuthHandler {
 public:
  explicit TestServerAuthHandler(const std::string& username,
                                 const std::string& password);
  ~TestServerAuthHandler() override;
  Status Authenticate(const ServerCallContext& context, ServerAuthSender* outgoing,
                      ServerAuthReader* incoming) override;
  Status IsValid(const ServerCallContext& context, const std::string& token,
                 std::string* peer_identity) override;

 private:
  std::string username_;
  std::string password_;
};

class ARROW_FLIGHT_EXPORT TestServerBasicAuthHandler : public ServerAuthHandler {
 public:
  explicit TestServerBasicAuthHandler(const std::string& username,
                                      const std::string& password);
  ~TestServerBasicAuthHandler() override;
  Status Authenticate(const ServerCallContext& context, ServerAuthSender* outgoing,
                      ServerAuthReader* incoming) override;
  Status IsValid(const ServerCallContext& context, const std::string& token,
                 std::string* peer_identity) override;

 private:
  BasicAuth basic_auth_;
};

class ARROW_FLIGHT_EXPORT TestClientAuthHandler : public ClientAuthHandler {
 public:
  explicit TestClientAuthHandler(const std::string& username,
                                 const std::string& password);
  ~TestClientAuthHandler() override;
  Status Authenticate(ClientAuthSender* outgoing, ClientAuthReader* incoming) override;
  Status GetToken(std::string* token) override;

 private:
  std::string username_;
  std::string password_;
};

class ARROW_FLIGHT_EXPORT TestClientBasicAuthHandler : public ClientAuthHandler {
 public:
  explicit TestClientBasicAuthHandler(const std::string& username,
                                      const std::string& password);
  ~TestClientBasicAuthHandler() override;
  Status Authenticate(ClientAuthSender* outgoing, ClientAuthReader* incoming) override;
  Status GetToken(std::string* token) override;

 private:
  BasicAuth basic_auth_;
  std::string token_;
};

}  // namespace arrow::flight
