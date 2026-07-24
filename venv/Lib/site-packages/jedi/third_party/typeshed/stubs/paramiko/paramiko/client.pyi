from _typeshed import FileDescriptorOrPath
from collections.abc import Iterable, Mapping
from typing import NoReturn, Protocol, type_check_only

from paramiko.auth_strategy import AuthStrategy
from paramiko.channel import Channel, ChannelFile, ChannelStderrFile, ChannelStdinFile
from paramiko.hostkeys import HostKeys
from paramiko.pkey import PKey
from paramiko.sftp_client import SFTPClient
from paramiko.transport import Transport, _SocketLike
from paramiko.util import ClosingContextManager

@type_check_only
class _TransportFactory(Protocol):
    def __call__(
        self,
        sock: _SocketLike,
        /,
        *,
        gss_kex: bool,
        gss_deleg_creds: bool,
        disabled_algorithms: Mapping[str, Iterable[str]] | None,
    ) -> Transport: ...

class SSHClient(ClosingContextManager):
    def __init__(self) -> None: ...
    def load_system_host_keys(self, filename: FileDescriptorOrPath | None = None) -> None: ...
    def load_host_keys(self, filename: FileDescriptorOrPath) -> None: ...
    def save_host_keys(self, filename: FileDescriptorOrPath) -> None: ...
    def get_host_keys(self) -> HostKeys: ...
    def set_log_channel(self, name: str) -> None: ...
    def set_missing_host_key_policy(self, policy: type[MissingHostKeyPolicy] | MissingHostKeyPolicy) -> None: ...
    def connect(
        self,
        hostname: str,
        port: int = 22,
        username: str | None = None,
        password: str | None = None,
        pkey: PKey | None = None,
        key_filename: str | None = None,
        timeout: float | None = None,
        allow_agent: bool = True,
        look_for_keys: bool = True,
        compress: bool = False,
        sock: _SocketLike | None = None,
        gss_auth: bool = False,
        gss_kex: bool = False,
        gss_deleg_creds: bool = True,
        gss_host: str | None = None,
        banner_timeout: float | None = None,
        auth_timeout: float | None = None,
        channel_timeout: float | None = None,
        gss_trust_dns: bool = True,
        passphrase: str | None = None,
        disabled_algorithms: Mapping[str, Iterable[str]] | None = None,
        transport_factory: _TransportFactory | None = None,
        auth_strategy: AuthStrategy | None = None,
    ) -> None: ...
    def close(self) -> None: ...
    def exec_command(
        self,
        command: str,
        bufsize: int = -1,
        timeout: float | None = None,
        get_pty: bool = False,
        environment: Mapping[str, str] | None = None,
    ) -> tuple[ChannelStdinFile, ChannelFile, ChannelStderrFile]: ...
    def invoke_shell(
        self,
        term: str = "vt100",
        width: int = 80,
        height: int = 24,
        width_pixels: int = 0,
        height_pixels: int = 0,
        environment: Mapping[str, str] | None = None,
    ) -> Channel: ...
    def open_sftp(self) -> SFTPClient: ...
    def get_transport(self) -> Transport | None: ...

class MissingHostKeyPolicy:
    def missing_host_key(self, client: SSHClient, hostname: str, key: PKey) -> None: ...

class AutoAddPolicy(MissingHostKeyPolicy):
    def missing_host_key(self, client: SSHClient, hostname: str, key: PKey) -> None: ...

class RejectPolicy(MissingHostKeyPolicy):
    def missing_host_key(self, client: SSHClient, hostname: str, key: PKey) -> NoReturn: ...

class WarningPolicy(MissingHostKeyPolicy):
    def missing_host_key(self, client: SSHClient, hostname: str, key: PKey) -> None: ...
