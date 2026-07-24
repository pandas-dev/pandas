from typing import Any

class XRaySession: ...
class XRayQuery: ...

class XRaySessionMaker:
    def __init__(
        self,
        bind=None,
        class_: type = ...,
        autoflush: bool = True,
        autocommit: bool = False,
        expire_on_commit: bool = True,
        info: dict[Any, Any] | None = None,  # it was taken from sqlalchemy stubs
        **kw,
    ) -> None: ...
