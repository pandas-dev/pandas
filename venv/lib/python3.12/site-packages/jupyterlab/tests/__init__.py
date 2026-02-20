# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from typing import NamedTuple


class Response(NamedTuple):
    """Fake tornado response."""

    body: bytes


def fake_client_factory():
    class FakeClient:
        """Fake AsyncHTTPClient

        body can be set in the test to a custom value.
        """

        body = b""

        async def fetch(*args, **kwargs):
            return Response(FakeClient.body)

    return FakeClient
