from _typeshed import Incomplete

from gunicorn.glogging import Logger as GLogger

from .._types import _ASGIAppType

class LifespanManager:
    app: _ASGIAppType
    logger: GLogger
    state: dict[Incomplete, Incomplete]

    def __init__(self, app: _ASGIAppType, logger: GLogger, state: dict[Incomplete, Incomplete] | None = None) -> None: ...
    async def startup(self) -> None: ...
    async def shutdown(self) -> None: ...
