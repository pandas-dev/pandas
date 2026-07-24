from _typeshed import Unused
from _typeshed.wsgi import WSGIApplication
from collections.abc import Callable, Iterable, Mapping
from typing import Any, Literal

from waitress.adjustments import _AdjustmentsParams
from waitress.server import BaseWSGIServer

def serve(
    app: WSGIApplication,
    *,
    _server: Callable[..., BaseWSGIServer] = ...,
    _quiet: bool = False,
    _profile: bool = False,
    **kw: _AdjustmentsParams,
) -> None: ...
def serve_paste(app: WSGIApplication, global_conf: Unused, **kw: _AdjustmentsParams) -> Literal[0]: ...
def profile(cmd: str, globals: dict[str, Any], locals: Mapping[str, Any], sort_order: Iterable[str], callers: bool) -> None: ...
