from _typeshed import Incomplete
from typing import ClassVar

class Environment:
    Development: ClassVar[Environment]
    QA: ClassVar[Environment]
    Sandbox: ClassVar[Environment]
    Production: ClassVar[Environment]
    All: ClassVar[dict[str, Environment]]
    __name__: str
    is_ssl: bool
    ssl_certificate: Incomplete
    def __init__(
        self,
        name,
        server: str,
        port,
        auth_url: str,
        is_ssl: bool,
        ssl_certificate,
        graphql_server: str = "",
        graphql_port: str = "",
    ) -> None: ...
    @property
    def base_url(self) -> str: ...
    @property
    def port(self) -> int: ...
    @property
    def auth_url(self) -> str: ...
    @property
    def protocol(self) -> str: ...
    @property
    def server(self) -> str: ...
    @property
    def server_and_port(self) -> str: ...
    @property
    def graphql_server(self) -> str: ...
    @property
    def graphql_port(self) -> str: ...
    @property
    def graphql_server_and_port(self) -> str: ...
    @staticmethod
    def parse_environment(environment: Environment | str | None) -> Environment | None: ...
    @staticmethod
    def braintree_root() -> str: ...
