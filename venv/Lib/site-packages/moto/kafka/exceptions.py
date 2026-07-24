"""Exceptions raised by the kafka service."""

from moto.core.exceptions import JsonRESTError


class KafkaException(JsonRESTError):
    pass


class NotFoundException(KafkaException):
    def __init__(self, msg: str):
        super().__init__("NotFoundException", msg)


class BadRequestException(KafkaException):
    def __init__(self, msg: str):
        super().__init__("BadRequestException", msg)
