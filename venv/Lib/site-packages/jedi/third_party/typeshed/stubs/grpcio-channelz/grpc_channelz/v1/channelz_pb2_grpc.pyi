from _typeshed import Incomplete
from typing import Final

import grpc

GRPC_GENERATED_VERSION: Final[str]
GRPC_VERSION: Final[str]

class ChannelzStub:
    GetTopChannels: Incomplete
    GetServers: Incomplete
    GetServer: Incomplete
    GetServerSockets: Incomplete
    GetChannel: Incomplete
    GetSubchannel: Incomplete
    GetSocket: Incomplete
    def __init__(self, channel: grpc.Channel): ...

class ChannelzServicer:
    def GetTopChannels(self, request, context): ...
    def GetServers(self, request, context): ...
    def GetServer(self, request, context): ...
    def GetServerSockets(self, request, context): ...
    def GetChannel(self, request, context): ...
    def GetSubchannel(self, request, context): ...
    def GetSocket(self, request, context): ...

def add_ChannelzServicer_to_server(servicer, server) -> None: ...

class Channelz:
    @staticmethod
    def GetTopChannels(
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
    def GetServers(
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
    def GetServer(
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
    def GetServerSockets(
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
    def GetChannel(
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
    def GetSubchannel(
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
    def GetSocket(
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
