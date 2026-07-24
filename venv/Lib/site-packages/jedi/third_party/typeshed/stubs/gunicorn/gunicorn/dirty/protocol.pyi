import asyncio
import socket
from _typeshed import Incomplete
from typing import ClassVar, Final

MAGIC: Final = b"GD"
VERSION: Final = 0x01
MSG_TYPE_REQUEST: Final = 0x01
MSG_TYPE_RESPONSE: Final = 0x02
MSG_TYPE_ERROR: Final = 0x03
MSG_TYPE_CHUNK: Final = 0x04
MSG_TYPE_END: Final = 0x05
MSG_TYPE_STASH: Final = 0x10
MSG_TYPE_STATUS: Final = 0x11
MSG_TYPE_MANAGE: Final = 0x12
MSG_TYPE_REQUEST_STR: Final = "request"
MSG_TYPE_RESPONSE_STR: Final = "response"
MSG_TYPE_ERROR_STR: Final = "error"
MSG_TYPE_CHUNK_STR: Final = "chunk"
MSG_TYPE_END_STR: Final = "end"
MSG_TYPE_STASH_STR: Final = "stash"
MSG_TYPE_STATUS_STR: Final = "status"
MSG_TYPE_MANAGE_STR: Final = "manage"
MSG_TYPE_TO_STR: Final[dict[int, str]]
MSG_TYPE_FROM_STR: Final[dict[str, int]]
STASH_OP_PUT: Final = 1
STASH_OP_GET: Final = 2
STASH_OP_DELETE: Final = 3
STASH_OP_KEYS: Final = 4
STASH_OP_CLEAR: Final = 5
STASH_OP_INFO: Final = 6
STASH_OP_ENSURE: Final = 7
STASH_OP_DELETE_TABLE: Final = 8
STASH_OP_TABLES: Final = 9
STASH_OP_EXISTS: Final = 10
MANAGE_OP_ADD: Final = 1
MANAGE_OP_REMOVE: Final = 2
HEADER_FORMAT: Final = ">2sBBIQ"
HEADER_SIZE: Final[int]
MAX_MESSAGE_SIZE: Final = 67108864

class BinaryProtocol:
    HEADER_SIZE: ClassVar[int]
    MAX_MESSAGE_SIZE: ClassVar[int]
    MSG_TYPE_REQUEST: ClassVar[str]
    MSG_TYPE_RESPONSE: ClassVar[str]
    MSG_TYPE_ERROR: ClassVar[str]
    MSG_TYPE_CHUNK: ClassVar[str]
    MSG_TYPE_END: ClassVar[str]
    MSG_TYPE_STASH: ClassVar[str]
    MSG_TYPE_STATUS: ClassVar[str]
    MSG_TYPE_MANAGE: ClassVar[str]

    @staticmethod
    def encode_header(msg_type: int, request_id: int, payload_length: int) -> bytes: ...
    @staticmethod
    def decode_header(data: bytes) -> tuple[int, int, int]: ...
    @staticmethod
    def encode_request(
        request_id: int,
        app_path: str,
        action: str,
        args: tuple[Incomplete, ...] | None = None,
        kwargs: dict[str, Incomplete] | None = None,
    ) -> bytes: ...
    @staticmethod
    def encode_response(request_id: int, result) -> bytes: ...
    @staticmethod
    def encode_error(request_id: int, error: BaseException | dict[str, Incomplete]) -> bytes: ...
    @staticmethod
    def encode_chunk(request_id: int, data) -> bytes: ...
    @staticmethod
    def encode_end(request_id: int) -> bytes: ...
    @staticmethod
    def encode_status(request_id: int) -> bytes: ...
    @staticmethod
    def encode_manage(request_id: int, op: int, count: int = 1) -> bytes: ...
    @staticmethod
    def encode_stash(request_id: int, op: int, table: str, key=None, value=None, pattern=None) -> bytes: ...
    @staticmethod
    def decode_message(data: bytes) -> tuple[str, int, Incomplete]: ...
    @staticmethod
    async def read_message_async(reader: asyncio.StreamReader) -> dict[str, Incomplete]: ...
    @staticmethod
    async def write_message_async(writer: asyncio.StreamWriter, message: dict[str, Incomplete]) -> None: ...
    @staticmethod
    def _recv_exactly(sock: socket.socket, n: int) -> bytes: ...
    @staticmethod
    def read_message(sock: socket.socket) -> dict[str, Incomplete]: ...
    @staticmethod
    def write_message(sock: socket.socket, message: dict[str, Incomplete]) -> None: ...
    @staticmethod
    def _encode_from_dict(message: dict[str, Incomplete]) -> bytes: ...

DirtyProtocol = BinaryProtocol

# TODO: Use TypedDict for results
def make_request(
    request_id: int | str,
    app_path: str,
    action: str,
    args: tuple[Incomplete, ...] | None = None,
    kwargs: dict[str, Incomplete] | None = None,
) -> dict[str, Incomplete]: ...
def make_response(request_id: int | str, result) -> dict[str, Incomplete]: ...
def make_error_response(request_id: int | str, error) -> dict[str, Incomplete]: ...
def make_chunk_message(request_id: int | str, data) -> dict[str, Incomplete]: ...
def make_end_message(request_id: int | str) -> dict[str, Incomplete]: ...
def make_stash_message(
    request_id: int | str, op: int, table: str, key=None, value=None, pattern=None
) -> dict[str, Incomplete]: ...
def make_manage_message(request_id: int | str, op: int, count: int = 1) -> dict[str, Incomplete]: ...
