from typing import Any

from moto.core.exceptions import ServiceException


class STSClientError(ServiceException):
    pass


class STSValidationError(STSClientError):
    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__("ValidationError", *args, **kwargs)
