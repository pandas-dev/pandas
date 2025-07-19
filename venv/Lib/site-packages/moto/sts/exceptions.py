from typing import Any

from moto.core.exceptions import RESTError


class STSClientError(RESTError):
    code = 400


class STSValidationError(STSClientError):
    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__("ValidationError", *args, **kwargs)
