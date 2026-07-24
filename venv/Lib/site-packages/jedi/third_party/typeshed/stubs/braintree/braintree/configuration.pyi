from _typeshed import Incomplete

from braintree.braintree_gateway import BraintreeGateway
from braintree.util.graphql_client import GraphQLClient
from braintree.util.http import Http

class Configuration:
    @staticmethod
    def configure(
        environment,
        merchant_id: str,
        public_key: str,
        private_key: str,
        *,
        http_strategy=None,
        timeout: int = 60,
        wrap_http_exceptions: bool = False,
    ) -> None: ...
    @staticmethod
    def for_partner(
        environment,
        partner_id: str,
        public_key: str,
        private_key: str,
        *,
        http_strategy=None,
        timeout: int = 60,
        wrap_http_exceptions: bool = False,
    ) -> Configuration: ...
    @staticmethod
    def gateway() -> BraintreeGateway: ...
    @staticmethod
    def instantiate() -> Configuration: ...
    @staticmethod
    def api_version() -> str: ...
    @staticmethod
    def graphql_api_version() -> str: ...
    environment: Incomplete
    merchant_id: str | None
    public_key: str | None
    private_key: str | None
    client_id: str | None
    client_secret: str | None
    access_token: str | None
    timeout: int
    wrap_http_exceptions: bool
    def __init__(
        self,
        environment=None,
        merchant_id: str | None = None,
        public_key: str | None = None,
        private_key: str | None = None,
        client_id: str | None = None,
        client_secret: str | None = None,
        access_token: str | None = None,
        *args,
        timeout: int = 60,
        wrap_http_exceptions: bool = False,
        http_strategy=None,
    ) -> None: ...
    def base_merchant_path(self) -> str: ...
    def base_url(self) -> str: ...
    def graphql_base_url(self) -> str: ...
    def http(self) -> Http: ...
    def graphql_client(self) -> GraphQLClient: ...
    def http_strategy(self): ...
    def close(self) -> None: ...
    def has_client_credentials(self) -> bool: ...
    def assert_has_client_credentials(self) -> None: ...
    def has_access_token(self) -> bool: ...
