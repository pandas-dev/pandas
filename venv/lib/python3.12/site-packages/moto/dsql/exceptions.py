"""Exceptions raised by the dsql service."""

from moto.core.exceptions import JsonRESTError


class ValidationException(JsonRESTError):
    """Tag validation failed."""

    code = 400

    def __init__(self, message: str):
        super().__init__("ValidationException", message)
