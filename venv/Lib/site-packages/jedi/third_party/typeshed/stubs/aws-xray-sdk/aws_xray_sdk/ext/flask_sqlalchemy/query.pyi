from _typeshed import Incomplete

from aws_xray_sdk.ext.sqlalchemy.query import XRaySession

class XRayBaseQuery: ...

class XRaySignallingSession(XRaySession):
    app: Incomplete
    def __init__(self, db, autocommit: bool = False, autoflush: bool = True, **options) -> None: ...
    def get_bind(self, mapper=None, clause=None): ...

class XRayFlaskSqlAlchemy:
    def __init__(
        self,
        app=None,
        use_native_unicode: bool = True,
        session_options=None,
        metadata=None,
        query_class: type = ...,
        model_class: type = ...,
    ) -> None: ...
    def create_session(self, options: dict[str, Incomplete]): ...
