from _typeshed import Incomplete
from collections.abc import Iterable
from typing import Final
from typing_extensions import TypeAlias

import grpc
import grpc.aio
from google.protobuf import descriptor_pool
from grpc_reflection.v1alpha import reflection_pb2 as _reflection_pb2
from grpc_reflection.v1alpha._base import BaseReflectionServicer

from . import _async as aio

SERVICE_NAME: Final[str]

_AnyServer: TypeAlias = grpc.Server | grpc.aio.Server
_AnyServicerContext: TypeAlias = grpc.ServicerContext | grpc.aio.ServicerContext[Incomplete, Incomplete]

class ReflectionServicer(BaseReflectionServicer):
    def ServerReflectionInfo(
        self, request_iterator: Iterable[_reflection_pb2.ServerReflectionRequest], context: _AnyServicerContext
    ): ...

def enable_server_reflection(
    service_names: Iterable[str], server: _AnyServer, pool: descriptor_pool.DescriptorPool | None = None
) -> None: ...

__all__ = ["SERVICE_NAME", "ReflectionServicer", "enable_server_reflection", "aio"]
