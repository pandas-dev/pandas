from _collections_abc import Generator, dict_keys
from _typeshed import Incomplete, ReadableBuffer
from types import TracebackType
from typing import Literal
from typing_extensions import Self, TypeAlias

from pyasn1.type.base import Asn1Item

from .pooling import ServerPool
from .server import Server

SASL_AVAILABLE_MECHANISMS: Incomplete
CLIENT_STRATEGIES: Incomplete

_ServerSequence: TypeAlias = (
    set[Server] | list[Server] | tuple[Server, ...] | Generator[Server, None, None] | dict_keys[Server, Incomplete]
)

class Connection:
    connection_lock: Incomplete
    last_error: str
    strategy_type: Incomplete
    user: Incomplete
    password: Incomplete
    authentication: Incomplete
    version: Incomplete
    auto_referrals: Incomplete
    request: Incomplete
    response: Incomplete | None
    result: Incomplete
    bound: bool
    listening: bool
    closed: bool
    auto_bind: Incomplete
    sasl_mechanism: Incomplete
    sasl_credentials: Incomplete
    socket: Incomplete
    tls_started: bool
    sasl_in_progress: bool
    read_only: Incomplete
    lazy: Incomplete
    pool_name: Incomplete
    pool_size: int | None
    cred_store: Incomplete
    pool_lifetime: Incomplete
    pool_keepalive: Incomplete
    starting_tls: bool
    check_names: Incomplete
    raise_exceptions: Incomplete
    auto_range: Incomplete
    extend: Incomplete
    fast_decoder: Incomplete
    receive_timeout: Incomplete
    empty_attributes: Incomplete
    use_referral_cache: Incomplete
    auto_escape: Incomplete
    auto_encode: Incomplete
    source_address: Incomplete
    source_port_list: Incomplete
    server_pool: Incomplete | None
    server: Incomplete
    strategy: Incomplete
    send: Incomplete
    open: Incomplete
    get_response: Incomplete
    post_send_single_response: Incomplete
    post_send_search: Incomplete
    def __init__(
        self,
        server: Server | str | _ServerSequence | ServerPool,
        user: str | None = None,
        password: str | None = None,
        auto_bind: Literal["DEFAULT", "NONE", "NO_TLS", "TLS_BEFORE_BIND", "TLS_AFTER_BIND"] | bool = "DEFAULT",
        version: int = 3,
        authentication: Literal["ANONYMOUS", "SIMPLE", "SASL", "NTLM"] | None = None,
        client_strategy: Literal[
            "SYNC",
            "SAFE_RESTARTABLE",
            "SAFE_SYNC",
            "ASYNC",
            "LDIF",
            "RESTARTABLE",
            "REUSABLE",
            "MOCK_SYNC",
            "MOCK_ASYNC",
            "ASYNC_STREAM",
        ] = "SYNC",
        auto_referrals: bool = True,
        auto_range: bool = True,
        sasl_mechanism: str | None = None,
        sasl_credentials=None,
        check_names: bool = True,
        collect_usage: bool = False,
        read_only: bool = False,
        lazy: bool = False,
        raise_exceptions: bool = False,
        pool_name: str | None = None,
        pool_size: int | None = None,
        pool_lifetime: int | None = None,
        cred_store=None,
        fast_decoder: bool = True,
        receive_timeout=None,
        return_empty_attributes: bool = True,
        use_referral_cache: bool = False,
        auto_escape: bool = True,
        auto_encode: bool = True,
        pool_keepalive=None,
        source_address: str | None = None,
        source_port: int | None = None,
        source_port_list=None,
    ) -> None: ...
    def repr_with_sensitive_data_stripped(self): ...
    @property
    def stream(self): ...
    @stream.setter
    def stream(self, value) -> None: ...
    @property
    def usage(self): ...
    def __enter__(self) -> Self: ...
    def __exit__(
        self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> Literal[False] | None: ...
    def bind(self, read_server_info: bool = True, controls=None): ...
    def rebind(
        self,
        user=None,
        password=None,
        authentication=None,
        sasl_mechanism=None,
        sasl_credentials=None,
        read_server_info: bool = True,
        controls=None,
    ): ...
    def unbind(self, controls=None): ...
    def search(
        self,
        search_base: str,
        search_filter: str,
        search_scope: Literal["BASE", "LEVEL", "SUBTREE"] = "SUBTREE",
        dereference_aliases: Literal["NEVER", "SEARCH", "FINDING_BASE", "ALWAYS"] = "ALWAYS",
        attributes=None,
        size_limit: int = 0,
        time_limit: int = 0,
        types_only: bool = False,
        get_operational_attributes: bool = False,
        controls=None,
        paged_size: int | None = None,
        paged_criticality: bool = False,
        paged_cookie: str | bytes | None = None,
        auto_escape: bool | None = None,
    ): ...
    def compare(self, dn, attribute, value, controls=None): ...
    def add(self, dn, object_class=None, attributes=None, controls=None): ...
    def delete(self, dn, controls=None): ...
    def modify(self, dn, changes, controls=None): ...
    def modify_dn(self, dn, relative_dn, delete_old_dn: bool = True, new_superior=None, controls=None): ...
    def abandon(self, message_id, controls=None): ...
    def extended(
        self, request_name, request_value: Asn1Item | ReadableBuffer | None = None, controls=None, no_encode: bool | None = None
    ): ...
    def start_tls(self, read_server_info: bool = True): ...
    def do_sasl_bind(self, controls): ...
    def do_ntlm_bind(self, controls): ...
    def refresh_server_info(self) -> None: ...
    def response_to_ldif(
        self, search_result=None, all_base64: bool = False, line_separator=None, sort_order=None, stream=None
    ): ...
    def response_to_json(
        self,
        raw: bool = False,
        search_result=None,
        indent: int = 4,
        sort: bool = True,
        stream=None,
        checked_attributes: bool = True,
        include_empty: bool = True,
    ): ...
    def response_to_file(self, target, raw: bool = False, indent: int = 4, sort: bool = True) -> None: ...
    @property
    def entries(self): ...
