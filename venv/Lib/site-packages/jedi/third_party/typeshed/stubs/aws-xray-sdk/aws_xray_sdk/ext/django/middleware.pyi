from _typeshed import Incomplete
from logging import Logger
from typing import Final

log: Logger

USER_AGENT_KEY: Final = "HTTP_USER_AGENT"
X_FORWARDED_KEY: Final = "HTTP_X_FORWARDED_FOR"
REMOTE_ADDR_KEY: Final = "REMOTE_ADDR"
HOST_KEY: Final = "HTTP_HOST"
CONTENT_LENGTH_KEY: Final = "content-length"

class XRayMiddleware:
    get_response: Incomplete
    in_lambda_ctx: bool

    def __init__(self, get_response) -> None: ...
    def __call__(self, request): ...
    def process_exception(self, request, exception: Exception) -> None: ...
