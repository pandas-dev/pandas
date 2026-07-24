"""Exceptions raised by the timestreaminfluxdb service."""

from moto.core.exceptions import JsonRESTError


class TimestreamInfluxDBException(JsonRESTError):
    pass


class ValidationException(TimestreamInfluxDBException):
    code = 400

    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class ConflictException(TimestreamInfluxDBException):
    code = 400

    def __init__(self, message: str):
        super().__init__("ConflictException", message)


class ResourceNotFoundException(TimestreamInfluxDBException):
    code = 400

    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)
