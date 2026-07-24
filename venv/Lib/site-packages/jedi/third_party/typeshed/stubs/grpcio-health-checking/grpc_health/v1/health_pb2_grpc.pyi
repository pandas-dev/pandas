from typing import Final

GRPC_GENERATED_VERSION: Final[str]
GRPC_VERSION: Final[str]

class HealthStub:
    def __init__(self, channel) -> None: ...

class HealthServicer:
    def Check(self, request, context): ...
    def Watch(self, request, context): ...

def add_HealthServicer_to_server(servicer, server) -> None: ...

class Health:
    @staticmethod
    def Check(
        request,
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
    @staticmethod
    def Watch(
        request,
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
