from grpc_channelz.v1 import _async as aio
from grpc_channelz.v1._servicer import ChannelzServicer

def add_channelz_servicer(server) -> None: ...

__all__ = ["aio", "add_channelz_servicer", "ChannelzServicer"]
