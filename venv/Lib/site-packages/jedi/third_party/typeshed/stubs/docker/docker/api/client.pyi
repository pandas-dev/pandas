from _typeshed import Incomplete
from collections.abc import Mapping, Sequence

import requests
from docker.tls import TLSConfig
from requests.adapters import BaseAdapter

from .build import BuildApiMixin
from .config import ConfigApiMixin
from .container import ContainerApiMixin
from .daemon import DaemonApiMixin
from .exec_api import ExecApiMixin
from .image import ImageApiMixin
from .network import NetworkApiMixin
from .plugin import PluginApiMixin
from .secret import SecretApiMixin
from .service import ServiceApiMixin
from .swarm import SwarmApiMixin
from .volume import VolumeApiMixin

class APIClient(
    requests.Session,
    BuildApiMixin,
    ConfigApiMixin,
    ContainerApiMixin,
    DaemonApiMixin,
    ExecApiMixin,
    ImageApiMixin,
    NetworkApiMixin,
    PluginApiMixin,
    SecretApiMixin,
    ServiceApiMixin,
    SwarmApiMixin,
    VolumeApiMixin,
):
    __attrs__: Sequence[str]
    base_url: str
    timeout: int
    credstore_env: Mapping[Incomplete, Incomplete] | None
    def __init__(
        self,
        base_url: str | None = None,
        version: str | None = None,
        timeout: int = 60,
        tls: bool | TLSConfig = False,
        user_agent: str = ...,
        num_pools: int | None = None,
        credstore_env: Mapping[Incomplete, Incomplete] | None = None,
        use_ssh_client: bool = False,
        max_pool_size: int = 10,
    ) -> None: ...
    def get_adapter(self, url: str) -> BaseAdapter: ...
    @property
    def api_version(self) -> str: ...
    def reload_config(self, dockercfg_path: str | None = None) -> None: ...
