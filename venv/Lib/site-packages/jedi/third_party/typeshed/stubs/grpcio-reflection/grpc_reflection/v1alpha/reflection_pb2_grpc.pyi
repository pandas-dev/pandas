from binascii import Incomplete
from typing import Final

import grpc

GRPC_GENERATED_VERSION: Final[str]
GRPC_VERSION: Final[str]

class ServerReflectionStub:
    ServerReflectionInfo: Incomplete
    def __init__(self, channel: grpc.Channel) -> None: ...

class ServerReflectionServicer:
    def ServerReflectionInfo(self, request_iterator, context): ...

def add_ServerReflectionServicer_to_server(servicer, server): ...

class ServerReflection:
    @staticmethod
    def ServerReflectionInfo(
        request_iterator,
        target,
        options=(),
        channel_credentials=None,
        call_credentials=None,
        insecure=False,
        compression=None,
        wait_for_ready=None,
        timeout=None,
        metadata=None,
    ): ...
