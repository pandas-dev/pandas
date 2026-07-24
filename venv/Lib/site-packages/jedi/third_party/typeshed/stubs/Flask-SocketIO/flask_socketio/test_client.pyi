from _typeshed import Incomplete
from typing import Any, TypedDict, type_check_only

from flask import Flask
from flask.testing import FlaskClient

@type_check_only
class _Packet(TypedDict):
    name: str
    args: Any
    namespace: str

class SocketIOTestClient:
    def __init__(
        self,
        app: Flask,
        socketio,
        namespace: str | None = None,
        query_string: str | None = None,
        headers: dict[str, Incomplete] | None = None,
        auth: dict[str, Incomplete] | None = None,
        flask_test_client: FlaskClient | None = None,
    ) -> None: ...
    def is_connected(self, namespace: str | None = None) -> bool: ...
    def connect(
        self,
        namespace: str | None = None,
        query_string: str | None = None,
        headers: dict[str, Incomplete] | None = None,
        auth: dict[str, Incomplete] | None = None,
    ) -> None: ...
    def disconnect(self, namespace: str | None = None) -> None: ...
    def emit(self, event: str, *args, callback: bool = False, namespace: str | None = None) -> Incomplete | None: ...
    def send(
        self,
        data: str | dict[str, Incomplete] | list[Incomplete],
        json: bool = False,
        callback: bool = False,
        namespace: str | None = None,
    ): ...
    def get_received(self, namespace: str | None = None) -> list[_Packet]: ...
