from _typeshed import Incomplete
from typing import Final

REQUESTS_SESSION_KWARGS: Final = ["proxies", "hooks", "stream", "verify", "cert", "max_redirects", "trust_env"]

def update_session_configure(session, kwargs: dict[str, Incomplete]) -> None: ...
