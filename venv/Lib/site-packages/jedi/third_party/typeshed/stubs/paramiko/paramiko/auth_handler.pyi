from collections.abc import Callable
from threading import Event
from typing_extensions import TypeAlias

from paramiko.message import Message
from paramiko.pkey import PKey
from paramiko.ssh_gss import _SSH_GSSAuth
from paramiko.transport import Transport

_InteractiveCallback: TypeAlias = Callable[[str, str, list[tuple[str, bool]]], list[str]]

class AuthHandler:
    transport: Transport
    username: str | None
    authenticated: bool
    auth_event: Event | None
    auth_method: str
    banner: str | None
    password: str | None
    private_key: PKey | None
    interactive_handler: _InteractiveCallback | None
    submethods: str | None
    auth_username: str | None
    auth_fail_count: int
    gss_host: str | None
    gss_deleg_creds: bool
    def __init__(self, transport: Transport) -> None: ...
    def is_authenticated(self) -> bool: ...
    def get_username(self) -> str | None: ...
    def auth_none(self, username: str, event: Event) -> None: ...
    def auth_publickey(self, username: str, key: PKey, event: Event) -> None: ...
    def auth_password(self, username: str, password: str, event: Event) -> None: ...
    def auth_interactive(self, username: str, handler: _InteractiveCallback, event: Event, submethods: str = "") -> None: ...
    def auth_gssapi_with_mic(self, username: str, gss_host: str, gss_deleg_creds: bool, event: Event) -> None: ...
    def auth_gssapi_keyex(self, username: str, event: Event) -> None: ...
    def abort(self) -> None: ...
    def wait_for_response(self, event: Event) -> list[str]: ...

class GssapiWithMicAuthHandler:
    method: str
    sshgss: _SSH_GSSAuth
    def __init__(self, delegate: AuthHandler, sshgss: _SSH_GSSAuth) -> None: ...
    def abort(self) -> None: ...
    @property
    def transport(self) -> Transport: ...
    @property
    def auth_username(self) -> str: ...
    @property
    def gss_host(self) -> str: ...

class AuthOnlyHandler(AuthHandler):
    def send_auth_request(
        self, username: str, method: str, finish_message: Callable[[Message], None] | None = None
    ) -> list[str]: ...
    def auth_none(self, username: str) -> list[str]: ...  # type: ignore[override]
    def auth_publickey(self, username: str, key: PKey) -> list[str]: ...  # type: ignore[override]
    def auth_password(self, username: str, password: str) -> list[str]: ...  # type: ignore[override]
    def auth_interactive(self, username: str, handler: _InteractiveCallback, submethods: str = "") -> list[str]: ...  # type: ignore[override]
