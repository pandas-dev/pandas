"""Exceptions raised by the cloudhsmv2 service."""

from moto.core.exceptions import JsonRESTError


class CloudHSMv2ClientError(JsonRESTError):
    """Base class for CloudHSMv2 errors."""

    code = 400


class ResourceNotFoundException(CloudHSMv2ClientError):
    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)


class InvalidRequestException(CloudHSMv2ClientError):
    def __init__(self, message: str):
        super().__init__("InvalidRequestException", message)


class ClientError(CloudHSMv2ClientError):
    def __init__(self, error_type: str, message: str):
        super().__init__(error_type, message)
