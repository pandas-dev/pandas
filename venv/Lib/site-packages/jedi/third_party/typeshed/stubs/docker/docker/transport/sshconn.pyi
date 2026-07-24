import socket
import subprocess

import urllib3
import urllib3.connection
from docker.transport.basehttpadapter import BaseHTTPAdapter
from paramiko import SSHClient, Transport
from urllib3._collections import RecentlyUsedContainer as urllib3_RecentlyUsedContainer

RecentlyUsedContainer = urllib3_RecentlyUsedContainer

class SSHSocket(socket.socket):
    host: str
    port: str | None
    user: str | None
    proc: subprocess.Popen[bytes] | None
    def __init__(self, host: str) -> None: ...
    def connect(self, **kwargs) -> None: ...  # type: ignore[override]
    def sendall(self, data) -> None: ...  # type: ignore[override]
    def send(self, data): ...  # type: ignore[override]
    def recv(self, n): ...  # type: ignore[override]
    def makefile(self, mode): ...  # type: ignore[override]
    def close(self) -> None: ...

class SSHConnection(urllib3.connection.HTTPConnection):
    ssh_transport: Transport | None
    timeout: int
    ssh_host: str | None
    def __init__(self, ssh_transport: Transport | None = None, timeout: int = 60, host: str | None = None) -> None: ...
    sock: SSHSocket | None
    def connect(self) -> None: ...

class SSHConnectionPool(urllib3.connectionpool.HTTPConnectionPool):
    scheme: str
    ssh_transport: Transport | None
    timeout: urllib3.Timeout
    ssh_host: str | None
    def __init__(
        self, ssh_client: SSHClient | None = None, timeout: int = 60, maxsize: int = 10, host: str | None = None
    ) -> None: ...

class SSHHTTPAdapter(BaseHTTPAdapter):
    __attrs__: list[str]
    ssh_client: SSHClient | None
    ssh_host: str
    timeout: int
    max_pool_size: int
    pools: int
    def __init__(
        self, base_url: str, timeout: int = 60, pool_connections: int = 25, max_pool_size: int = 10, shell_out: bool = False
    ) -> None: ...
    def get_connection(self, url: str | bytes, proxies=None) -> SSHConnectionPool: ...
    def close(self) -> None: ...
