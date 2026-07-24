"""BedrockAgentCoreControl exceptions."""

from moto.core.exceptions import ServiceException


class BedrockAgentCoreControlClientError(ServiceException):
    pass


class ResourceNotFoundException(BedrockAgentCoreControlClientError):
    code = "ResourceNotFoundException"


class ConflictException(BedrockAgentCoreControlClientError):
    code = "ConflictException"
