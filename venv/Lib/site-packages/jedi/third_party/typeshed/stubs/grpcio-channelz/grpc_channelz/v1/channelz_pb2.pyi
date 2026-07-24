from _typeshed import Incomplete
from collections.abc import Iterable, Mapping
from typing import ClassVar, final

from google._upb._message import Descriptor, FileDescriptor, MessageMeta
from google.protobuf import any_pb2, duration_pb2, message, timestamp_pb2, wrappers_pb2
from google.protobuf.internal import containers

DESCRIPTOR: FileDescriptor

@final
class Channel(message.Message, metaclass=MessageMeta):
    REF_FIELD_NUMBER: ClassVar[int]
    DATA_FIELD_NUMBER: ClassVar[int]
    CHANNEL_REF_FIELD_NUMBER: ClassVar[int]
    SUBCHANNEL_REF_FIELD_NUMBER: ClassVar[int]
    SOCKET_REF_FIELD_NUMBER: ClassVar[int]
    ref: ChannelRef
    data: ChannelData
    channel_ref: containers.RepeatedCompositeFieldContainer[ChannelRef]
    subchannel_ref: containers.RepeatedCompositeFieldContainer[SubchannelRef]
    socket_ref: containers.RepeatedCompositeFieldContainer[SocketRef]
    def __init__(
        self,
        ref: ChannelRef | Mapping[Incomplete, Incomplete] | None = ...,
        data: ChannelData | Mapping[Incomplete, Incomplete] | None = ...,
        channel_ref: Iterable[ChannelRef | Mapping[Incomplete, Incomplete]] | None = ...,
        subchannel_ref: Iterable[SubchannelRef | Mapping[Incomplete, Incomplete]] | None = ...,
        socket_ref: Iterable[SocketRef | Mapping[Incomplete, Incomplete]] | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class Subchannel(message.Message, metaclass=MessageMeta):
    REF_FIELD_NUMBER: ClassVar[int]
    DATA_FIELD_NUMBER: ClassVar[int]
    CHANNEL_REF_FIELD_NUMBER: ClassVar[int]
    SUBCHANNEL_REF_FIELD_NUMBER: ClassVar[int]
    SOCKET_REF_FIELD_NUMBER: ClassVar[int]
    ref: SubchannelRef
    data: ChannelData
    channel_ref: containers.RepeatedCompositeFieldContainer[ChannelRef]
    subchannel_ref: containers.RepeatedCompositeFieldContainer[SubchannelRef]
    socket_ref: containers.RepeatedCompositeFieldContainer[SocketRef]
    def __init__(
        self,
        ref: SubchannelRef | Mapping[Incomplete, Incomplete] | None = ...,
        data: ChannelData | Mapping[Incomplete, Incomplete] | None = ...,
        channel_ref: Iterable[ChannelRef | Mapping[Incomplete, Incomplete]] | None = ...,
        subchannel_ref: Iterable[SubchannelRef | Mapping[Incomplete, Incomplete]] | None = ...,
        socket_ref: Iterable[SocketRef | Mapping[Incomplete, Incomplete]] | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class ChannelConnectivityState(message.Message, metaclass=MessageMeta):
    State: Incomplete
    UNKNOWN: Incomplete
    IDLE: Incomplete
    CONNECTING: Incomplete
    READY: Incomplete
    TRANSIENT_FAILURE: Incomplete
    SHUTDOWN: Incomplete
    STATE_FIELD_NUMBER: ClassVar[int]
    state: Incomplete
    def __init__(self, state: Incomplete | str | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class ChannelData(message.Message, metaclass=MessageMeta):
    STATE_FIELD_NUMBER: ClassVar[int]
    TARGET_FIELD_NUMBER: ClassVar[int]
    TRACE_FIELD_NUMBER: ClassVar[int]
    CALLS_STARTED_FIELD_NUMBER: ClassVar[int]
    CALLS_SUCCEEDED_FIELD_NUMBER: ClassVar[int]
    CALLS_FAILED_FIELD_NUMBER: ClassVar[int]
    LAST_CALL_STARTED_TIMESTAMP_FIELD_NUMBER: ClassVar[int]
    state: ChannelConnectivityState
    target: str
    trace: ChannelTrace
    calls_started: int
    calls_succeeded: int
    calls_failed: int
    last_call_started_timestamp: timestamp_pb2.Timestamp
    def __init__(
        self,
        state: ChannelConnectivityState | Mapping[Incomplete, Incomplete] | None = ...,
        target: str | None = ...,
        trace: ChannelTrace | Mapping[Incomplete, Incomplete] | None = ...,
        calls_started: int | None = ...,
        calls_succeeded: int | None = ...,
        calls_failed: int | None = ...,
        last_call_started_timestamp: timestamp_pb2.Timestamp | Mapping[Incomplete, Incomplete] | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class ChannelTraceEvent(message.Message, metaclass=MessageMeta):
    Severity: Incomplete
    CT_UNKNOWN: Incomplete
    CT_INFO: Incomplete
    CT_WARNING: Incomplete
    CT_ERROR: Incomplete
    DESCRIPTION_FIELD_NUMBER: ClassVar[int]
    SEVERITY_FIELD_NUMBER: ClassVar[int]
    TIMESTAMP_FIELD_NUMBER: ClassVar[int]
    CHANNEL_REF_FIELD_NUMBER: ClassVar[int]
    SUBCHANNEL_REF_FIELD_NUMBER: ClassVar[int]
    description: str
    severity: Incomplete
    timestamp: timestamp_pb2.Timestamp
    channel_ref: ChannelRef
    subchannel_ref: SubchannelRef
    def __init__(
        self,
        description: str | None = ...,
        severity: Incomplete | str | None = ...,
        timestamp: timestamp_pb2.Timestamp | Mapping[Incomplete, Incomplete] | None = ...,
        channel_ref: ChannelRef | Mapping[Incomplete, Incomplete] | None = ...,
        subchannel_ref: SubchannelRef | Mapping[Incomplete, Incomplete] | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class ChannelTrace(message.Message, metaclass=MessageMeta):
    NUM_EVENTS_LOGGED_FIELD_NUMBER: ClassVar[int]
    CREATION_TIMESTAMP_FIELD_NUMBER: ClassVar[int]
    EVENTS_FIELD_NUMBER: ClassVar[int]
    num_events_logged: int
    creation_timestamp: timestamp_pb2.Timestamp
    events: containers.RepeatedCompositeFieldContainer[ChannelTraceEvent]
    def __init__(
        self,
        num_events_logged: int | None = ...,
        creation_timestamp: timestamp_pb2.Timestamp | Mapping[Incomplete, Incomplete] | None = ...,
        events: Iterable[ChannelTraceEvent | Mapping[Incomplete, Incomplete]] | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class ChannelRef(message.Message, metaclass=MessageMeta):
    CHANNEL_ID_FIELD_NUMBER: ClassVar[int]
    NAME_FIELD_NUMBER: ClassVar[int]
    channel_id: int
    name: str
    def __init__(self, channel_id: int | None = ..., name: str | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class SubchannelRef(message.Message, metaclass=MessageMeta):
    SUBCHANNEL_ID_FIELD_NUMBER: ClassVar[int]
    NAME_FIELD_NUMBER: ClassVar[int]
    subchannel_id: int
    name: str
    def __init__(self, subchannel_id: int | None = ..., name: str | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class SocketRef(message.Message, metaclass=MessageMeta):
    SOCKET_ID_FIELD_NUMBER: ClassVar[int]
    NAME_FIELD_NUMBER: ClassVar[int]
    socket_id: int
    name: str
    def __init__(self, socket_id: int | None = ..., name: str | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class ServerRef(message.Message, metaclass=MessageMeta):
    SERVER_ID_FIELD_NUMBER: ClassVar[int]
    NAME_FIELD_NUMBER: ClassVar[int]
    server_id: int
    name: str
    def __init__(self, server_id: int | None = ..., name: str | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class Server(message.Message, metaclass=MessageMeta):
    REF_FIELD_NUMBER: ClassVar[int]
    DATA_FIELD_NUMBER: ClassVar[int]
    LISTEN_SOCKET_FIELD_NUMBER: ClassVar[int]
    ref: ServerRef
    data: ServerData
    listen_socket: containers.RepeatedCompositeFieldContainer[SocketRef]
    def __init__(
        self,
        ref: ServerRef | Mapping[Incomplete, Incomplete] | None = ...,
        data: ServerData | Mapping[Incomplete, Incomplete] | None = ...,
        listen_socket: Iterable[SocketRef | Mapping[Incomplete, Incomplete]] | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class ServerData(message.Message, metaclass=MessageMeta):
    TRACE_FIELD_NUMBER: ClassVar[int]
    CALLS_STARTED_FIELD_NUMBER: ClassVar[int]
    CALLS_SUCCEEDED_FIELD_NUMBER: ClassVar[int]
    CALLS_FAILED_FIELD_NUMBER: ClassVar[int]
    LAST_CALL_STARTED_TIMESTAMP_FIELD_NUMBER: ClassVar[int]
    trace: ChannelTrace
    calls_started: int
    calls_succeeded: int
    calls_failed: int
    last_call_started_timestamp: timestamp_pb2.Timestamp
    def __init__(
        self,
        trace: ChannelTrace | Mapping[Incomplete, Incomplete] | None = ...,
        calls_started: int | None = ...,
        calls_succeeded: int | None = ...,
        calls_failed: int | None = ...,
        last_call_started_timestamp: timestamp_pb2.Timestamp | Mapping[Incomplete, Incomplete] | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class Socket(message.Message, metaclass=MessageMeta):
    REF_FIELD_NUMBER: ClassVar[int]
    DATA_FIELD_NUMBER: ClassVar[int]
    LOCAL_FIELD_NUMBER: ClassVar[int]
    REMOTE_FIELD_NUMBER: ClassVar[int]
    SECURITY_FIELD_NUMBER: ClassVar[int]
    REMOTE_NAME_FIELD_NUMBER: ClassVar[int]
    ref: SocketRef
    data: SocketData
    local: Address
    remote: Address
    security: Security
    remote_name: str
    def __init__(
        self,
        ref: SocketRef | Mapping[Incomplete, Incomplete] | None = ...,
        data: SocketData | Mapping[Incomplete, Incomplete] | None = ...,
        local: Address | Mapping[Incomplete, Incomplete] | None = ...,
        remote: Address | Mapping[Incomplete, Incomplete] | None = ...,
        security: Security | Mapping[Incomplete, Incomplete] | None = ...,
        remote_name: str | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class SocketData(message.Message, metaclass=MessageMeta):
    STREAMS_STARTED_FIELD_NUMBER: ClassVar[int]
    STREAMS_SUCCEEDED_FIELD_NUMBER: ClassVar[int]
    STREAMS_FAILED_FIELD_NUMBER: ClassVar[int]
    MESSAGES_SENT_FIELD_NUMBER: ClassVar[int]
    MESSAGES_RECEIVED_FIELD_NUMBER: ClassVar[int]
    KEEP_ALIVES_SENT_FIELD_NUMBER: ClassVar[int]
    LAST_LOCAL_STREAM_CREATED_TIMESTAMP_FIELD_NUMBER: ClassVar[int]
    LAST_REMOTE_STREAM_CREATED_TIMESTAMP_FIELD_NUMBER: ClassVar[int]
    LAST_MESSAGE_SENT_TIMESTAMP_FIELD_NUMBER: ClassVar[int]
    LAST_MESSAGE_RECEIVED_TIMESTAMP_FIELD_NUMBER: ClassVar[int]
    LOCAL_FLOW_CONTROL_WINDOW_FIELD_NUMBER: ClassVar[int]
    REMOTE_FLOW_CONTROL_WINDOW_FIELD_NUMBER: ClassVar[int]
    OPTION_FIELD_NUMBER: ClassVar[int]
    streams_started: int
    streams_succeeded: int
    streams_failed: int
    messages_sent: int
    messages_received: int
    keep_alives_sent: int
    last_local_stream_created_timestamp: timestamp_pb2.Timestamp
    last_remote_stream_created_timestamp: timestamp_pb2.Timestamp
    last_message_sent_timestamp: timestamp_pb2.Timestamp
    last_message_received_timestamp: timestamp_pb2.Timestamp
    local_flow_control_window: wrappers_pb2.Int64Value
    remote_flow_control_window: wrappers_pb2.Int64Value
    option: containers.RepeatedCompositeFieldContainer[SocketOption]
    def __init__(
        self,
        streams_started: int | None = ...,
        streams_succeeded: int | None = ...,
        streams_failed: int | None = ...,
        messages_sent: int | None = ...,
        messages_received: int | None = ...,
        keep_alives_sent: int | None = ...,
        last_local_stream_created_timestamp: timestamp_pb2.Timestamp | Mapping[Incomplete, Incomplete] | None = ...,
        last_remote_stream_created_timestamp: timestamp_pb2.Timestamp | Mapping[Incomplete, Incomplete] | None = ...,
        last_message_sent_timestamp: timestamp_pb2.Timestamp | Mapping[Incomplete, Incomplete] | None = ...,
        last_message_received_timestamp: timestamp_pb2.Timestamp | Mapping[Incomplete, Incomplete] | None = ...,
        local_flow_control_window: wrappers_pb2.Int64Value | Mapping[Incomplete, Incomplete] | None = ...,
        remote_flow_control_window: wrappers_pb2.Int64Value | Mapping[Incomplete, Incomplete] | None = ...,
        option: Iterable[SocketOption | Mapping[Incomplete, Incomplete]] | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class Address(message.Message, metaclass=MessageMeta):
    @final
    class TcpIpAddress(message.Message, metaclass=MessageMeta):
        IP_ADDRESS_FIELD_NUMBER: ClassVar[int]
        PORT_FIELD_NUMBER: ClassVar[int]
        ip_address: bytes
        port: int
        def __init__(self, ip_address: bytes | None = ..., port: int | None = ...) -> None: ...

    @final
    class UdsAddress(message.Message, metaclass=MessageMeta):
        FILENAME_FIELD_NUMBER: ClassVar[int]
        filename: str
        def __init__(self, filename: str | None = ...) -> None: ...

    @final
    class OtherAddress(message.Message, metaclass=MessageMeta):
        NAME_FIELD_NUMBER: ClassVar[int]
        VALUE_FIELD_NUMBER: ClassVar[int]
        name: str
        value: any_pb2.Any
        def __init__(self, name: str | None = ..., value: any_pb2.Any | Mapping[Incomplete, Incomplete] | None = ...) -> None: ...

    TCPIP_ADDRESS_FIELD_NUMBER: ClassVar[int]
    UDS_ADDRESS_FIELD_NUMBER: ClassVar[int]
    OTHER_ADDRESS_FIELD_NUMBER: ClassVar[int]
    tcpip_address: Address.TcpIpAddress
    uds_address: Address.UdsAddress
    other_address: Address.OtherAddress
    def __init__(
        self,
        tcpip_address: Address.TcpIpAddress | Mapping[Incomplete, Incomplete] | None = ...,
        uds_address: Address.UdsAddress | Mapping[Incomplete, Incomplete] | None = ...,
        other_address: Address.OtherAddress | Mapping[Incomplete, Incomplete] | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class Security(message.Message, metaclass=MessageMeta):
    @final
    class Tls(message.Message, metaclass=MessageMeta):
        STANDARD_NAME_FIELD_NUMBER: ClassVar[int]
        OTHER_NAME_FIELD_NUMBER: ClassVar[int]
        LOCAL_CERTIFICATE_FIELD_NUMBER: ClassVar[int]
        REMOTE_CERTIFICATE_FIELD_NUMBER: ClassVar[int]
        standard_name: str
        other_name: str
        local_certificate: bytes
        remote_certificate: bytes
        def __init__(
            self,
            standard_name: str | None = ...,
            other_name: str | None = ...,
            local_certificate: bytes | None = ...,
            remote_certificate: bytes | None = ...,
        ) -> None: ...

    @final
    class OtherSecurity(message.Message, metaclass=MessageMeta):
        NAME_FIELD_NUMBER: ClassVar[int]
        VALUE_FIELD_NUMBER: ClassVar[int]
        name: str
        value: any_pb2.Any
        def __init__(self, name: str | None = ..., value: any_pb2.Any | Mapping[Incomplete, Incomplete] | None = ...) -> None: ...

    TLS_FIELD_NUMBER: ClassVar[int]
    OTHER_FIELD_NUMBER: ClassVar[int]
    tls: Security.Tls
    other: Security.OtherSecurity
    def __init__(
        self,
        tls: Security.Tls | Mapping[Incomplete, Incomplete] | None = ...,
        other: Security.OtherSecurity | Mapping[Incomplete, Incomplete] | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class SocketOption(message.Message, metaclass=MessageMeta):
    NAME_FIELD_NUMBER: ClassVar[int]
    VALUE_FIELD_NUMBER: ClassVar[int]
    ADDITIONAL_FIELD_NUMBER: ClassVar[int]
    name: str
    value: str
    additional: any_pb2.Any
    def __init__(
        self,
        name: str | None = ...,
        value: str | None = ...,
        additional: any_pb2.Any | Mapping[Incomplete, Incomplete] | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class SocketOptionTimeout(message.Message, metaclass=MessageMeta):
    DURATION_FIELD_NUMBER: ClassVar[int]
    duration: duration_pb2.Duration
    def __init__(self, duration: duration_pb2.Duration | Mapping[Incomplete, Incomplete] | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class SocketOptionLinger(message.Message, metaclass=MessageMeta):
    ACTIVE_FIELD_NUMBER: ClassVar[int]
    DURATION_FIELD_NUMBER: ClassVar[int]
    active: bool
    duration: duration_pb2.Duration
    def __init__(
        self, active: bool = ..., duration: duration_pb2.Duration | Mapping[Incomplete, Incomplete] | None = ...
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class SocketOptionTcpInfo(message.Message, metaclass=MessageMeta):
    TCPI_STATE_FIELD_NUMBER: ClassVar[int]
    TCPI_CA_STATE_FIELD_NUMBER: ClassVar[int]
    TCPI_RETRANSMITS_FIELD_NUMBER: ClassVar[int]
    TCPI_PROBES_FIELD_NUMBER: ClassVar[int]
    TCPI_BACKOFF_FIELD_NUMBER: ClassVar[int]
    TCPI_OPTIONS_FIELD_NUMBER: ClassVar[int]
    TCPI_SND_WSCALE_FIELD_NUMBER: ClassVar[int]
    TCPI_RCV_WSCALE_FIELD_NUMBER: ClassVar[int]
    TCPI_RTO_FIELD_NUMBER: ClassVar[int]
    TCPI_ATO_FIELD_NUMBER: ClassVar[int]
    TCPI_SND_MSS_FIELD_NUMBER: ClassVar[int]
    TCPI_RCV_MSS_FIELD_NUMBER: ClassVar[int]
    TCPI_UNACKED_FIELD_NUMBER: ClassVar[int]
    TCPI_SACKED_FIELD_NUMBER: ClassVar[int]
    TCPI_LOST_FIELD_NUMBER: ClassVar[int]
    TCPI_RETRANS_FIELD_NUMBER: ClassVar[int]
    TCPI_FACKETS_FIELD_NUMBER: ClassVar[int]
    TCPI_LAST_DATA_SENT_FIELD_NUMBER: ClassVar[int]
    TCPI_LAST_ACK_SENT_FIELD_NUMBER: ClassVar[int]
    TCPI_LAST_DATA_RECV_FIELD_NUMBER: ClassVar[int]
    TCPI_LAST_ACK_RECV_FIELD_NUMBER: ClassVar[int]
    TCPI_PMTU_FIELD_NUMBER: ClassVar[int]
    TCPI_RCV_SSTHRESH_FIELD_NUMBER: ClassVar[int]
    TCPI_RTT_FIELD_NUMBER: ClassVar[int]
    TCPI_RTTVAR_FIELD_NUMBER: ClassVar[int]
    TCPI_SND_SSTHRESH_FIELD_NUMBER: ClassVar[int]
    TCPI_SND_CWND_FIELD_NUMBER: ClassVar[int]
    TCPI_ADVMSS_FIELD_NUMBER: ClassVar[int]
    TCPI_REORDERING_FIELD_NUMBER: ClassVar[int]
    tcpi_state: int
    tcpi_ca_state: int
    tcpi_retransmits: int
    tcpi_probes: int
    tcpi_backoff: int
    tcpi_options: int
    tcpi_snd_wscale: int
    tcpi_rcv_wscale: int
    tcpi_rto: int
    tcpi_ato: int
    tcpi_snd_mss: int
    tcpi_rcv_mss: int
    tcpi_unacked: int
    tcpi_sacked: int
    tcpi_lost: int
    tcpi_retrans: int
    tcpi_fackets: int
    tcpi_last_data_sent: int
    tcpi_last_ack_sent: int
    tcpi_last_data_recv: int
    tcpi_last_ack_recv: int
    tcpi_pmtu: int
    tcpi_rcv_ssthresh: int
    tcpi_rtt: int
    tcpi_rttvar: int
    tcpi_snd_ssthresh: int
    tcpi_snd_cwnd: int
    tcpi_advmss: int
    tcpi_reordering: int
    def __init__(
        self,
        tcpi_state: int | None = ...,
        tcpi_ca_state: int | None = ...,
        tcpi_retransmits: int | None = ...,
        tcpi_probes: int | None = ...,
        tcpi_backoff: int | None = ...,
        tcpi_options: int | None = ...,
        tcpi_snd_wscale: int | None = ...,
        tcpi_rcv_wscale: int | None = ...,
        tcpi_rto: int | None = ...,
        tcpi_ato: int | None = ...,
        tcpi_snd_mss: int | None = ...,
        tcpi_rcv_mss: int | None = ...,
        tcpi_unacked: int | None = ...,
        tcpi_sacked: int | None = ...,
        tcpi_lost: int | None = ...,
        tcpi_retrans: int | None = ...,
        tcpi_fackets: int | None = ...,
        tcpi_last_data_sent: int | None = ...,
        tcpi_last_ack_sent: int | None = ...,
        tcpi_last_data_recv: int | None = ...,
        tcpi_last_ack_recv: int | None = ...,
        tcpi_pmtu: int | None = ...,
        tcpi_rcv_ssthresh: int | None = ...,
        tcpi_rtt: int | None = ...,
        tcpi_rttvar: int | None = ...,
        tcpi_snd_ssthresh: int | None = ...,
        tcpi_snd_cwnd: int | None = ...,
        tcpi_advmss: int | None = ...,
        tcpi_reordering: int | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetTopChannelsRequest(message.Message, metaclass=MessageMeta):
    START_CHANNEL_ID_FIELD_NUMBER: ClassVar[int]
    MAX_RESULTS_FIELD_NUMBER: ClassVar[int]
    start_channel_id: int
    max_results: int
    def __init__(self, start_channel_id: int | None = ..., max_results: int | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetTopChannelsResponse(message.Message, metaclass=MessageMeta):
    CHANNEL_FIELD_NUMBER: ClassVar[int]
    END_FIELD_NUMBER: ClassVar[int]
    channel: containers.RepeatedCompositeFieldContainer[Channel]
    end: bool
    def __init__(self, channel: Iterable[Channel | Mapping[Incomplete, Incomplete]] | None = ..., end: bool = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetServersRequest(message.Message, metaclass=MessageMeta):
    START_SERVER_ID_FIELD_NUMBER: ClassVar[int]
    MAX_RESULTS_FIELD_NUMBER: ClassVar[int]
    start_server_id: int
    max_results: int
    def __init__(self, start_server_id: int | None = ..., max_results: int | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetServersResponse(message.Message, metaclass=MessageMeta):
    SERVER_FIELD_NUMBER: ClassVar[int]
    END_FIELD_NUMBER: ClassVar[int]
    server: containers.RepeatedCompositeFieldContainer[Server]
    end: bool
    def __init__(self, server: Iterable[Server | Mapping[Incomplete, Incomplete]] | None = ..., end: bool = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetServerRequest(message.Message, metaclass=MessageMeta):
    SERVER_ID_FIELD_NUMBER: ClassVar[int]
    server_id: int
    def __init__(self, server_id: int | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetServerResponse(message.Message, metaclass=MessageMeta):
    SERVER_FIELD_NUMBER: ClassVar[int]
    server: Server
    def __init__(self, server: Server | Mapping[Incomplete, Incomplete] | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetServerSocketsRequest(message.Message, metaclass=MessageMeta):
    SERVER_ID_FIELD_NUMBER: ClassVar[int]
    START_SOCKET_ID_FIELD_NUMBER: ClassVar[int]
    MAX_RESULTS_FIELD_NUMBER: ClassVar[int]
    server_id: int
    start_socket_id: int
    max_results: int
    def __init__(self, server_id: int | None = ..., start_socket_id: int | None = ..., max_results: int | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetServerSocketsResponse(message.Message, metaclass=MessageMeta):
    SOCKET_REF_FIELD_NUMBER: ClassVar[int]
    END_FIELD_NUMBER: ClassVar[int]
    socket_ref: containers.RepeatedCompositeFieldContainer[SocketRef]
    end: bool
    def __init__(
        self, socket_ref: Iterable[SocketRef | Mapping[Incomplete, Incomplete]] | None = ..., end: bool = ...
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetChannelRequest(message.Message, metaclass=MessageMeta):
    CHANNEL_ID_FIELD_NUMBER: ClassVar[int]
    channel_id: int
    def __init__(self, channel_id: int | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetChannelResponse(message.Message, metaclass=MessageMeta):
    CHANNEL_FIELD_NUMBER: ClassVar[int]
    channel: Channel
    def __init__(self, channel: Channel | Mapping[Incomplete, Incomplete] | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetSubchannelRequest(message.Message, metaclass=MessageMeta):
    SUBCHANNEL_ID_FIELD_NUMBER: ClassVar[int]
    subchannel_id: int
    def __init__(self, subchannel_id: int | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetSubchannelResponse(message.Message, metaclass=MessageMeta):
    SUBCHANNEL_FIELD_NUMBER: ClassVar[int]
    subchannel: Subchannel
    def __init__(self, subchannel: Subchannel | Mapping[Incomplete, Incomplete] | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetSocketRequest(message.Message, metaclass=MessageMeta):
    SOCKET_ID_FIELD_NUMBER: ClassVar[int]
    SUMMARY_FIELD_NUMBER: ClassVar[int]
    socket_id: int
    summary: bool
    def __init__(self, socket_id: int | None = ..., summary: bool = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class GetSocketResponse(message.Message, metaclass=MessageMeta):
    SOCKET_FIELD_NUMBER: ClassVar[int]
    socket: Socket
    def __init__(self, socket: Socket | Mapping[Incomplete, Incomplete] | None = ...) -> None: ...
    DESCRIPTOR: Descriptor
