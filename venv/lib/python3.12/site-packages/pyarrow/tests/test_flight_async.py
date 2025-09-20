# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

import asyncio

import pytest

import pyarrow

flight = pytest.importorskip("pyarrow.flight")
pytestmark = pytest.mark.flight


class ExampleServer(flight.FlightServerBase):
    simple_info = flight.FlightInfo(
        pyarrow.schema([("a", "int32")]),
        flight.FlightDescriptor.for_command(b"simple"),
        []
    )

    def get_flight_info(self, context, descriptor):
        if descriptor.command == b"simple":
            return self.simple_info
        elif descriptor.command == b"unknown":
            raise NotImplementedError("Unknown command")

        raise NotImplementedError("Unknown descriptor")


def async_or_skip(client):
    if not client.supports_async:
        # Use async error message as skip message
        with pytest.raises(NotImplementedError) as e:
            client.as_async()
        pytest.skip(str(e.value))


@pytest.fixture(scope="module")
def flight_client():
    with ExampleServer() as server:
        with flight.connect(f"grpc://localhost:{server.port}") as client:
            yield client


@pytest.fixture(scope="module")
def async_client(flight_client):
    async_or_skip(flight_client)
    yield flight_client.as_async()


def test_async_support_property(flight_client):
    assert isinstance(flight_client.supports_async, bool)
    if flight_client.supports_async:
        flight_client.as_async()
    else:
        with pytest.raises(NotImplementedError):
            flight_client.as_async()


def test_get_flight_info(async_client):
    async def _test():
        descriptor = flight.FlightDescriptor.for_command(b"simple")
        info = await async_client.get_flight_info(descriptor)
        assert info == ExampleServer.simple_info

    asyncio.run(_test())


def test_get_flight_info_error(async_client):
    async def _test():
        descriptor = flight.FlightDescriptor.for_command(b"unknown")
        with pytest.raises(NotImplementedError) as excinfo:
            await async_client.get_flight_info(descriptor)

        assert "Unknown command" in repr(excinfo.value)

    asyncio.run(_test())
