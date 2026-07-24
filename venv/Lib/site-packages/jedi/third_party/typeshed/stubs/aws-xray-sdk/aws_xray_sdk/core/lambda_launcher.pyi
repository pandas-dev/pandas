from logging import Logger
from typing import Final

from .context import Context

log: Logger
LAMBDA_TRACE_HEADER_KEY: Final = "_X_AMZN_TRACE_ID"
LAMBDA_TASK_ROOT_KEY: Final = "LAMBDA_TASK_ROOT"
TOUCH_FILE_DIR: Final = "/tmp/.aws-xray/"
TOUCH_FILE_PATH: Final = "/tmp/.aws-xray/initialized"

def check_in_lambda() -> LambdaContext | None: ...

class LambdaContext(Context):
    def __init__(self) -> None: ...
    @property  # type: ignore[override]
    def context_missing(self) -> None: ...
    @context_missing.setter
    def context_missing(self, value: str) -> None: ...
