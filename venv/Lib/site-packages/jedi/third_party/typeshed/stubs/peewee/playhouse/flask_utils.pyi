from _typeshed import Unused
from collections.abc import Container
from typing import Any
from typing_extensions import TypeAlias

from peewee import Database, ModelBase, Proxy

# Is actually flask.Flask
_Flask: TypeAlias = Any

class FlaskDB:
    # Omitting undocumented base_model_class on purpose, use FlaskDB.Model instead
    database: Database | Proxy
    def __init__(
        self,
        app: _Flask | None = None,
        database: Database | Proxy | None = None,
        # Is actually type[ModelClass] but stubtest likely confuses with Model property
        # https://github.com/python/typeshed/pull/11731#issuecomment-2067694259
        model_class=...,
        excluded_routes: Container[str] | None = None,
    ) -> None: ...
    def init_app(self, app: _Flask) -> None: ...
    def get_model_class(self) -> type[ModelBase]: ...
    @property
    def Model(self) -> type[ModelBase]: ...
    def connect_db(self) -> None: ...
    def close_db(self, exc: Unused) -> None: ...
