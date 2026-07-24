import grpc_channelz.v1.channelz_pb2 as _channelz_pb2
import grpc_channelz.v1.channelz_pb2_grpc as _channelz_pb2_grpc
from grpc import ServicerContext

class ChannelzServicer(_channelz_pb2_grpc.ChannelzServicer):
    @staticmethod
    def GetTopChannels(
        request: _channelz_pb2.GetTopChannelsRequest, context: ServicerContext
    ) -> _channelz_pb2.GetTopChannelsResponse: ...
    @staticmethod
    def GetServers(request: _channelz_pb2.GetServersRequest, context: ServicerContext) -> _channelz_pb2.GetServersResponse: ...
    @staticmethod
    def GetServer(request: _channelz_pb2.GetServerRequest, context: ServicerContext) -> _channelz_pb2.GetServerResponse: ...
    @staticmethod
    def GetServerSockets(
        request: _channelz_pb2.GetServerSocketsRequest, context: ServicerContext
    ) -> _channelz_pb2.GetServerSocketsResponse: ...
    @staticmethod
    def GetChannel(request: _channelz_pb2.GetChannelRequest, context: ServicerContext) -> _channelz_pb2.GetChannelResponse: ...
    @staticmethod
    def GetSubchannel(
        request: _channelz_pb2.GetSubchannelRequest, context: ServicerContext
    ) -> _channelz_pb2.GetSubchannelResponse: ...
    @staticmethod
    def GetSocket(request: _channelz_pb2.GetSocketRequest, context: ServicerContext) -> _channelz_pb2.GetSocketResponse: ...
