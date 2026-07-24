import urllib3
import urllib3.connection
from docker.transport.basehttpadapter import BaseHTTPAdapter
from docker.transport.npipesocket import NpipeSocket
from urllib3._collections import RecentlyUsedContainer as urllib3_RecentlyUsedContainer

RecentlyUsedContainer = urllib3_RecentlyUsedContainer

class NpipeHTTPConnection(urllib3.connection.HTTPConnection):
    npipe_path: str
    timeout: int
    def __init__(self, npipe_path: str, timeout: int = 60) -> None: ...
    sock: NpipeSocket | None
    def connect(self) -> None: ...

class NpipeHTTPConnectionPool(urllib3.connectionpool.HTTPConnectionPool):
    npipe_path: str
    timeout: urllib3.Timeout
    def __init__(self, npipe_path: str, timeout: int = 60, maxsize: int = 10) -> None: ...

class NpipeHTTPAdapter(BaseHTTPAdapter):
    __attrs__: list[str]
    npipe_path: str
    timeout: int
    max_pool_size: int
    pools: RecentlyUsedContainer
    def __init__(self, base_url: str, timeout: int = 60, pool_connections: int = 25, max_pool_size: int = 10) -> None: ...
    def get_connection(self, url, proxies=None): ...
    def request_url(self, request, proxies): ...
