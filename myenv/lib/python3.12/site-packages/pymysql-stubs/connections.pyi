from _typeshed import Incomplete
from collections.abc import Mapping
from socket import socket as _socket
from ssl import _PasswordType
from typing import Any, AnyStr, Generic, TypeVar, overload
from typing_extensions import Self

from .charset import charset_by_id as charset_by_id, charset_by_name as charset_by_name
from .constants import CLIENT as CLIENT, COMMAND as COMMAND, FIELD_TYPE as FIELD_TYPE, SERVER_STATUS as SERVER_STATUS
from .cursors import Cursor
from .util import byte2int as byte2int, int2byte as int2byte

SSL_ENABLED: Any
DEFAULT_USER: Any
DEBUG: Any
DEFAULT_CHARSET: Any

_C = TypeVar("_C", bound=Cursor)
_C2 = TypeVar("_C2", bound=Cursor)

def dump_packet(data): ...
def pack_int24(n): ...
def lenenc_int(i: int) -> bytes: ...

class MysqlPacket:
    connection: Any
    def __init__(self, data, encoding): ...
    def get_all_data(self): ...
    def read(self, size): ...
    def read_all(self): ...
    def advance(self, length): ...
    def rewind(self, position: int = 0): ...
    def get_bytes(self, position, length: int = 1): ...
    def read_string(self) -> bytes: ...
    def read_uint8(self) -> Any: ...
    def read_uint16(self) -> Any: ...
    def read_uint24(self) -> Any: ...
    def read_uint32(self) -> Any: ...
    def read_uint64(self) -> Any: ...
    def read_length_encoded_integer(self) -> int: ...
    def read_length_coded_string(self) -> bytes: ...
    def read_struct(self, fmt: str) -> tuple[Any, ...]: ...
    def is_ok_packet(self) -> bool: ...
    def is_eof_packet(self) -> bool: ...
    def is_auth_switch_request(self) -> bool: ...
    def is_extra_auth_data(self) -> bool: ...
    def is_resultset_packet(self) -> bool: ...
    def is_load_local_packet(self) -> bool: ...
    def is_error_packet(self) -> bool: ...
    def check_error(self): ...
    def raise_for_error(self) -> None: ...
    def dump(self): ...

class FieldDescriptorPacket(MysqlPacket):
    def __init__(self, data, encoding): ...
    def description(self): ...
    def get_column_length(self): ...

class Connection(Generic[_C]):
    ssl: Any
    host: Any
    port: Any
    user: Any
    password: Any
    db: Any
    unix_socket: Any
    charset: str
    collation: str | None
    bind_address: Any
    use_unicode: Any
    client_flag: Any
    cursorclass: Any
    connect_timeout: Any
    messages: Any
    encoders: Any
    decoders: Any
    host_info: Any
    sql_mode: Any
    init_command: Any
    max_allowed_packet: int
    server_public_key: bytes
    @overload
    def __init__(
        self: Connection[Cursor],  # different between overloads
        *,
        host: str | None = None,
        user: Incomplete | None = None,
        password: str = "",
        database: Incomplete | None = None,
        port: int = 0,
        unix_socket: Incomplete | None = None,
        charset: str = "",
        collation: str | None = None,
        sql_mode: Incomplete | None = None,
        read_default_file: Incomplete | None = None,
        conv=None,
        use_unicode: bool | None = True,
        client_flag: int = 0,
        cursorclass: None = None,  # different between overloads
        init_command: Incomplete | None = None,
        connect_timeout: int | None = 10,
        ssl: Mapping[Any, Any] | None = None,
        ssl_ca=None,
        ssl_cert=None,
        ssl_disabled=None,
        ssl_key=None,
        ssl_key_password: _PasswordType | None = None,
        ssl_verify_cert=None,
        ssl_verify_identity=None,
        read_default_group: Incomplete | None = None,
        compress: Incomplete | None = None,
        named_pipe: Incomplete | None = None,
        autocommit: bool | None = False,
        db: Incomplete | None = None,
        passwd: Incomplete | None = None,
        local_infile: Incomplete | None = False,
        max_allowed_packet: int = 16777216,
        defer_connect: bool | None = False,
        auth_plugin_map: Mapping[Any, Any] | None = None,
        read_timeout: float | None = None,
        write_timeout: float | None = None,
        bind_address: Incomplete | None = None,
        binary_prefix: bool | None = False,
        program_name: Incomplete | None = None,
        server_public_key: bytes | None = None,
    ) -> None: ...
    @overload
    def __init__(
        # different between overloads:
        self: Connection[_C],  # pyright: ignore[reportInvalidTypeVarUse]  #11780
        *,
        host: str | None = None,
        user: Incomplete | None = None,
        password: str = "",
        database: Incomplete | None = None,
        port: int = 0,
        unix_socket: Incomplete | None = None,
        charset: str = "",
        collation: str | None = None,
        sql_mode: Incomplete | None = None,
        read_default_file: Incomplete | None = None,
        conv=None,
        use_unicode: bool | None = True,
        client_flag: int = 0,
        cursorclass: type[_C] = ...,  # different between overloads
        init_command: Incomplete | None = None,
        connect_timeout: int | None = 10,
        ssl: Mapping[Any, Any] | None = None,
        ssl_ca=None,
        ssl_cert=None,
        ssl_disabled=None,
        ssl_key=None,
        ssl_verify_cert=None,
        ssl_verify_identity=None,
        read_default_group: Incomplete | None = None,
        compress: Incomplete | None = None,
        named_pipe: Incomplete | None = None,
        autocommit: bool | None = False,
        db: Incomplete | None = None,
        passwd: Incomplete | None = None,
        local_infile: Incomplete | None = False,
        max_allowed_packet: int = 16777216,
        defer_connect: bool | None = False,
        auth_plugin_map: Mapping[Any, Any] | None = None,
        read_timeout: float | None = None,
        write_timeout: float | None = None,
        bind_address: Incomplete | None = None,
        binary_prefix: bool | None = False,
        program_name: Incomplete | None = None,
        server_public_key: bytes | None = None,
    ) -> None: ...
    socket: Any
    rfile: Any
    wfile: Any
    def close(self) -> None: ...
    @property
    def open(self) -> bool: ...
    def autocommit(self, value) -> None: ...
    def get_autocommit(self) -> bool: ...
    def commit(self) -> None: ...
    def begin(self) -> None: ...
    def rollback(self) -> None: ...
    def select_db(self, db) -> None: ...
    def escape(self, obj, mapping: Mapping[Any, Any] | None = None): ...
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
    def set_charset(self, charset) -> None: ...
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
    Warning: Any
    Error: Any
    InterfaceError: Any
    DatabaseError: Any
    DataError: Any
    OperationalError: Any
    IntegrityError: Any
    InternalError: Any
    ProgrammingError: Any
    NotSupportedError: Any

class MySQLResult:
    connection: Any
    affected_rows: Any
    insert_id: Any
    server_status: Any
    warning_count: Any
    message: Any
    field_count: Any
    description: Any
    rows: Any
    has_next: Any
    def __init__(self, connection: Connection[Any]) -> None: ...
    first_packet: Any
    def read(self) -> None: ...
    def init_unbuffered_query(self) -> None: ...

class LoadLocalFile:
    filename: Any
    connection: Connection[Any]
    def __init__(self, filename: Any, connection: Connection[Any]) -> None: ...
    def send_data(self) -> None: ...
