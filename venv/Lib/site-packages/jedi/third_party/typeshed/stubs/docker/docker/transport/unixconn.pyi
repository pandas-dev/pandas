import socket

import urllib3
import urllib3.connection
from docker.transport.basehttpadapter import BaseHTTPAdapter
from requests import PreparedRequest
from urllib3._collections import RecentlyUsedContainer as urllib3_RecentlyUsedContainer

RecentlyUsedContainer = urllib3_RecentlyUsedContainer

class UnixHTTPConnection(urllib3.connection.HTTPConnection):
    base_url: str
    unix_socket: str
    timeout: int
    def __init__(self, base_url: str, unix_socket: str, timeout: int = 60) -> None: ...
    sock: socket.socket | None
    def connect(self) -> None: ...

class UnixHTTPConnectionPool(urllib3.connectionpool.HTTPConnectionPool):
    base_url: str
    socket_path: str
    timeout: urllib3.Timeout
    def __init__(self, base_url: str, socket_path: str, timeout: int = 60, maxsize: int = 10) -> None: ...

class UnixHTTPAdapter(BaseHTTPAdapter):
    __attrs__: list[str]
    socket_path: str
    timeout: int
    max_pool_size: int
    pools: RecentlyUsedContainer
    def __init__(self, socket_url: str, timeout: int = 60, pool_connections: int = 25, max_pool_size: int = 10) -> None: ...
    def get_connection(self, url: bytes | str, proxies=None) -> UnixHTTPConnectionPool: ...
    # proxies is unused
    def request_url(self, request: PreparedRequest, proxies) -> str: ...
