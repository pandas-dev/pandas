"""Exceptions raised by the timestreamwrite service."""

from moto.core.exceptions import JsonRESTError


class ResourceNotFound(JsonRESTError):
    error_type = "com.amazonaws.timestream.v20181101#ResourceNotFoundException"

    def __init__(self, msg: str):
        super().__init__(ResourceNotFound.error_type, msg)
