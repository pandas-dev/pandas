from _typeshed import FileDescriptorOrPath, Incomplete, Unused
from collections.abc import Callable, Mapping
from socket import _Address, socket as _socket
from ssl import SSLContext, _PasswordType
from typing import Any, AnyStr, Generic, overload
from typing_extensions import Self, TypeVar, deprecated

from .charset import charset_by_id as charset_by_id, charset_by_name as charset_by_name
from .constants import CLIENT as CLIENT, COMMAND as COMMAND, FIELD_TYPE as FIELD_TYPE, SERVER_STATUS as SERVER_STATUS
from .cursors import Cursor
from .err import (
    DatabaseError,
    DataError,
    Error,
    IntegrityError,
    InterfaceError,
    InternalError,
    NotSupportedError,
    OperationalError,
    ProgrammingError,
    Warning,
)

_C = TypeVar("_C", bound=Cursor, default=Cursor)
_C2 = TypeVar("_C2", bound=Cursor)

SSL_ENABLED: bool
DEFAULT_USER: str | None
DEBUG: bool
DEFAULT_CHARSET: str
TEXT_TYPES: set[int]
MAX_PACKET_LEN: int

def dump_packet(data): ...
def _lenenc_int(i: int) -> bytes: ...

class Connection(Generic[_C]):
    ssl: bool
    host: str
    port: int
    user: str | bytes | None
    password: bytes
    db: str | bytes | None
    unix_socket: _Address | None
    charset: str
    collation: str | None
    bind_address: str | None
    use_unicode: bool
    client_flag: int
    cursorclass: type[_C]
    connect_timeout: float | None
    host_info: str
    sql_mode: str | None
    init_command: str | None
    max_allowed_packet: int
    server_public_key: bytes | None
    encoding: str
    autocommit_mode: bool | None
    encoders: dict[type[Any], Callable[[Any], str]]  # argument type depends on the key
    decoders: dict[int, Callable[[str], Any]]  # return type depends on the key

    @overload
    def __init__(
        self,
        *,
        user: str | bytes | None = None,
        password: str | bytes = "",
        host: str | None = None,
        database: str | bytes | None = None,
        unix_socket: _Address | None = None,
        port: int = 0,
        charset: str = "",
        collation: str | None = None,
        sql_mode: str | None = None,
        read_default_file: str | None = None,
        conv: dict[int | type[Any], Callable[[Any], str] | Callable[[str], Any]] | None = None,
        use_unicode: bool = True,
        client_flag: int = 0,
        cursorclass: type[_C] = ...,
        init_command: str | None = None,
        connect_timeout: float = 10,
        read_default_group: str | None = None,
        autocommit: bool | None = False,
        local_infile: bool = False,
        max_allowed_packet: int = 16_777_216,
        defer_connect: bool = False,
        auth_plugin_map: dict[str, Callable[[Connection[Any]], Any]] | None = None,
        read_timeout: float | None = None,
        write_timeout: float | None = None,
        bind_address: str | None = None,
        binary_prefix: bool = False,
        program_name: str | None = None,
        server_public_key: bytes | None = None,
        ssl: dict[str, Incomplete] | SSLContext | None = None,
        ssl_ca: str | None = None,
        ssl_cert: str | None = None,
        ssl_disabled: bool | None = None,
        ssl_key: str | None = None,
        ssl_key_password: _PasswordType | None = None,
        ssl_verify_cert: bool | None = None,
        ssl_verify_identity: bool | None = None,
        compress: Unused = None,
        named_pipe: Unused = None,
        # different between overloads:
        passwd: None = None,  # deprecated
        db: None = None,  # deprecated
    ) -> None: ...
    @overload
    @deprecated("'passwd' and 'db' arguments are deprecated. Use 'password' and 'database' instead.")
    def __init__(
        self,
        *,
        user: str | bytes | None = None,
        password: str | bytes = "",
        host: str | None = None,
        database: str | bytes | None = None,
        unix_socket: _Address | None = None,
        port: int = 0,
        charset: str = "",
        collation: str | None = None,
        sql_mode: str | None = None,
        read_default_file: str | None = None,
        conv: dict[int | type[Any], Callable[[Any], str] | Callable[[str], Any]] | None = None,
        use_unicode: bool = True,
        client_flag: int = 0,
        cursorclass: type[_C] = ...,
        init_command: str | None = None,
        connect_timeout: float = 10,
        read_default_group: str | None = None,
        autocommit: bool | None = False,
        local_infile: bool = False,
        max_allowed_packet: int = 16_777_216,
        defer_connect: bool = False,
        auth_plugin_map: dict[str, Callable[[Connection[Any]], Any]] | None = None,
        read_timeout: float | None = None,
        write_timeout: float | None = None,
        bind_address: str | None = None,
        binary_prefix: bool = False,
        program_name: str | None = None,
        server_public_key: bytes | None = None,
        ssl: dict[str, Incomplete] | SSLContext | None = None,
        ssl_ca: str | None = None,
        ssl_cert: str | None = None,
        ssl_disabled: bool | None = None,
        ssl_key: str | None = None,
        ssl_key_password: _PasswordType | None = None,
        ssl_verify_cert: bool | None = None,
        ssl_verify_identity: bool | None = None,
        compress: Unused = None,
        named_pipe: Unused = None,
        # different between overloads:
        passwd: str | bytes | None = None,  # deprecated
        db: str | bytes | None = None,  # deprecated
    ) -> None: ...
    def close(self) -> None: ...
    @property
    def open(self) -> bool: ...
    def __del__(self) -> None: ...
    def autocommit(self, value) -> None: ...
    def get_autocommit(self) -> bool: ...
    def commit(self) -> None: ...
    def begin(self) -> None: ...
    def rollback(self) -> None: ...
    def select_db(self, db) -> None: ...
    def escape(self, obj, mapping: Mapping[str, Incomplete] | None = None): ...
    def literal(self, obj): ...
    def escape_string(self, s: AnyStr) -> AnyStr: ...
    @overload
    def cursor(self, cursor: None = None) -> _C: ...
    @overload
    def cursor(self, cursor: type[_C2]) -> _C2: ...
    def query(self, sql, unbuffered: bool = False) -> int: ...
    def next_result(self, unbuffered: bool = False) -> int: ...
    def affected_rows(self): ...
    def kill(self, thread_id): ...
    def ping(self, reconnect: bool = True) -> None: ...
    @deprecated("Method is deprecated. Use set_character_set() instead.")
    def set_charset(self, charset: str) -> None: ...
    def set_character_set(self, charset: str, collation: str | None = None) -> None: ...
    def connect(self, sock: _socket | None = None) -> None: ...
    def write_packet(self, payload) -> None: ...
    def _read_packet(self, packet_type=...): ...
    def insert_id(self): ...
    def thread_id(self): ...
    def character_set_name(self): ...
    def get_host_info(self) -> str: ...
    def get_proto_info(self): ...
    def get_server_info(self): ...
    def show_warnings(self): ...
    def __enter__(self) -> Self: ...
    def __exit__(self, *exc_info: object) -> None: ...
    Warning: type[Warning]
    Error: type[Error]
    InterfaceError: type[InterfaceError]
    DatabaseError: type[DatabaseError]
    DataError: type[DataError]
    OperationalError: type[OperationalError]
    IntegrityError: type[IntegrityError]
    InternalError: type[InternalError]
    ProgrammingError: type[ProgrammingError]
    NotSupportedError: type[NotSupportedError]

class MySQLResult:
    connection: Connection[Any] | None
    affected_rows: int | None
    insert_id: int | None
    server_status: int | None
    warning_count: int
    message: str | None
    field_count: int
    description: Incomplete
    rows: Incomplete
    has_next: bool | None
    unbuffered_active: bool
    def __init__(self, connection: Connection[Any]) -> None: ...
    def __del__(self) -> None: ...
    first_packet: Incomplete
    def read(self) -> None: ...
    def init_unbuffered_query(self) -> None: ...

class LoadLocalFile:
    filename: FileDescriptorOrPath
    connection: Connection[Any]
    def __init__(self, filename: FileDescriptorOrPath, connection: Connection[Any]) -> None: ...
    def send_data(self) -> None: ...
