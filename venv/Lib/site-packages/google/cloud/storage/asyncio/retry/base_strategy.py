# Copyright 2025 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import abc
from typing import Any, Iterable


class _BaseResumptionStrategy(abc.ABC):
    """Abstract base class defining the interface for a bidi stream resumption strategy.

    This class defines the skeleton for a pluggable strategy that contains
    all the service-specific logic for a given bidi operation (e.g., reads
    or writes). This allows a generic retry manager to handle the common
    retry loop while sending the state management and request generation
    to a concrete implementation of this class.
    """

    @abc.abstractmethod
    def generate_requests(self, state: Any) -> Iterable[Any]:
        """Generates the next batch of requests based on the current state.

        This method is called at the beginning of each retry attempt. It should
        inspect the provided state object and generate the appropriate list of
        request protos to send to the server. For example, a read strategy
        would use this to implement "Smarter Resumption" by creating smaller
        `ReadRange` requests for partially downloaded ranges. For bidi-writes,
        it will set the `write_offset` field to the persisted size received
        from the server in the next request.

        :type state: Any
        :param state: An object containing all the state needed for the
                      operation (e.g., requested ranges, user buffers,
                      bytes written).
        """
        pass

    @abc.abstractmethod
    def update_state_from_response(self, response: Any, state: Any) -> None:
        """Updates the state based on a successful server response.

        This method is called for every message received from the server. It is
        responsible for processing the response and updating the shared state
        object.

        :type response: Any
        :param response: The response message received from the server.

        :type state: Any
        :param state: The shared state object for the operation, which will be
                      mutated by this method.
        """
        pass

    @abc.abstractmethod
    async def recover_state_on_failure(self, error: Exception, state: Any) -> None:
        """Prepares the state for the next retry attempt after a failure.

        This method is called when a retriable gRPC error occurs. It is
        responsible for performing any necessary actions to ensure the next
        retry attempt can succeed. For bidi reads, its primary role is to
        handle the `BidiReadObjectRedirectError` by extracting the
        `routing_token` and updating the state. For bidi writes, it will update
        the state to reflect any bytes that were successfully persisted before
        the failure.

        :type error: :class:`Exception`
        :param error: The exception that was caught by the retry engine.

        :type state: Any
        :param state: The shared state object for the operation.
        """
        pass
