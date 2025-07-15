from _typeshed import FileDescriptorOrPath, Incomplete
from collections.abc import Mapping
from socket import socket as _socket
from ssl import _PasswordType
from typing import AnyStr, Generic, TypeVar, overload
from typing_extensions import Self, deprecated

from .charset import charset_by_id as charset_by_id, charset_by_name as charset_by_name
from .constants import CLIENT as CLIENT, COMMAND as COMMAND, FIELD_TYPE as FIELD_TYPE, SERVER_STATUS as SERVER_STATUS
from .cursors import Cursor

_C = TypeVar("_C", bound=Cursor)
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
    ssl: Incomplete
    host: Incomplete
    port: Incomplete
    user: Incomplete
    password: Incomplete
    db: Incomplete
    unix_socket: Incomplete
    charset: str
    collation: str | None
    bind_address: Incomplete
    use_unicode: Incomplete
    client_flag: Incomplete
    cursorclass: Incomplete
    connect_timeout: Incomplete
    messages: Incomplete
    encoders: Incomplete
    decoders: Incomplete
    host_info: Incomplete
    sql_mode: Incomplete
    init_command: Incomplete
    max_allowed_packet: int
    server_public_key: bytes
    @overload
    def __init__(
        self: Connection[Cursor],  # different between overloads
        *,
        host: str | None = None,
        user=None,
        password: str = "",
        database=None,
        port: int = 0,
        unix_socket=None,
        charset: str = "",
        collation: str | None = None,
        sql_mode=None,
        read_default_file=None,
        conv=None,
        use_unicode: bool | None = True,
        client_flag: int = 0,
        cursorclass: None = None,  # different between overloads
        init_command=None,
        connect_timeout: int | None = 10,
        ssl: Mapping[Incomplete, Incomplete] | None = None,
        ssl_ca=None,
        ssl_cert=None,
        ssl_disabled=None,
        ssl_key=None,
        ssl_key_password: _PasswordType | None = None,
        ssl_verify_cert=None,
        ssl_verify_identity=None,
        read_default_group=None,
        compress=None,
        named_pipe=None,
        autocommit: bool | None = False,
        db=None,
        passwd=None,
        local_infile: Incomplete | None = False,
        max_allowed_packet: int = 16777216,
        defer_connect: bool | None = False,
        auth_plugin_map: Mapping[Incomplete, Incomplete] | None = None,
        read_timeout: float | None = None,
        write_timeout: float | None = None,
        bind_address=None,
        binary_prefix: bool | None = False,
        program_name=None,
        server_public_key: bytes | None = None,
    ) -> None: ...
    @overload
    def __init__(
        # different between overloads:
        self: Connection[_C],  # pyright: ignore[reportInvalidTypeVarUse]  #11780
        *,
        host: str | None = None,
        user=None,
        password: str = "",
        database=None,
        port: int = 0,
        unix_socket=None,
        charset: str = "",
        collation: str | None = None,
        sql_mode=None,
        read_default_file=None,
        conv=None,
        use_unicode: bool | None = True,
        client_flag: int = 0,
        cursorclass: type[_C] = ...,  # different between overloads
        init_command=None,
        connect_timeout: int | None = 10,
        ssl: Mapping[Incomplete, Incomplete] | None = None,
        ssl_ca=None,
        ssl_cert=None,
        ssl_disabled=None,
        ssl_key=None,
        ssl_verify_cert=None,
        ssl_verify_identity=None,
        read_default_group=None,
        compress=None,
        named_pipe=None,
        autocommit: bool | None = False,
        db=None,
        passwd=None,
        local_infile: Incomplete | None = False,
        max_allowed_packet: int = 16777216,
        defer_connect: bool | None = False,
        auth_plugin_map: Mapping[Incomplete, Incomplete] | None = None,
        read_timeout: float | None = None,
        write_timeout: float | None = None,
        bind_address=None,
        binary_prefix: bool | None = False,
        program_name=None,
        server_public_key: bytes | None = None,
    ) -> None: ...
    socket: Incomplete
    rfile: Incomplete
    wfile: Incomplete
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
    def escape(self, obj, mapping: Mapping[Incomplete, Incomplete] | None = None): ...
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
    def get_host_info(self): ...
    def get_proto_info(self): ...
    def get_server_info(self): ...
    def show_warnings(self): ...
    def __enter__(self) -> Self: ...
    def __exit__(self, *exc_info: object) -> None: ...
    Warning: Incomplete
    Error: Incomplete
    InterfaceError: Incomplete
    DatabaseError: Incomplete
    DataError: Incomplete
    OperationalError: Incomplete
    IntegrityError: Incomplete
    InternalError: Incomplete
    ProgrammingError: Incomplete
    NotSupportedError: Incomplete

class MySQLResult:
    connection: Incomplete
    affected_rows: Incomplete
    insert_id: Incomplete
    server_status: Incomplete
    warning_count: Incomplete
    message: Incomplete
    field_count: Incomplete
    description: Incomplete
    rows: Incomplete
    has_next: Incomplete
    def __init__(self, connection: Connection[Incomplete]) -> None: ...
    def __del__(self) -> None: ...
    first_packet: Incomplete
    def read(self) -> None: ...
    def init_unbuffered_query(self) -> None: ...

class LoadLocalFile:
    filename: FileDescriptorOrPath
    connection: Connection[Incomplete]
    def __init__(self, filename: FileDescriptorOrPath, connection: Connection[Incomplete]) -> None: ...
    def send_data(self) -> None: ...
