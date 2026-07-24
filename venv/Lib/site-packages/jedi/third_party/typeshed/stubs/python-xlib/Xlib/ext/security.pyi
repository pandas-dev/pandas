from _typeshed import Unused
from typing import Final

from Xlib.display import Display
from Xlib.protocol import rq
from Xlib.xobject import resource

extname: Final = "SECURITY"
SecurityClientTrusted: Final = 0
SecurityClientUntrusted: Final = 1
SecurityAuthorizationRevokedMask: Final = 1
AUTHID = rq.Card32

class QueryVersion(rq.ReplyRequest): ...

def query_version(self: Display | resource.Resource) -> QueryVersion: ...

class SecurityGenerateAuthorization(rq.ReplyRequest): ...

def generate_authorization(
    self: Display | resource.Resource,
    auth_proto: str,
    auth_data: bytes | bytearray = b"",
    timeout: int | None = None,
    trust_level: int | None = None,
    group: int | None = None,
    event_mask: int | None = None,
) -> SecurityGenerateAuthorization: ...

class SecurityRevokeAuthorization(rq.Request): ...

def revoke_authorization(self: Display | resource.Resource, authid: int) -> SecurityRevokeAuthorization: ...
def init(disp: Display, info: Unused) -> None: ...
