"""Exceptions raised by the bedrockruntime service."""

from moto.core.exceptions import ServiceException


class BedrockRuntimeError(ServiceException):
    pass
