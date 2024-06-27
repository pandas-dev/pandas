from moto.core.exceptions import AWSError


class InvalidRequestException(AWSError):
    TYPE = "InvalidRequestException"


class InvalidParameterValueException(AWSError):
    TYPE = "InvalidParameterValue"


class ValidationError(AWSError):
    TYPE = "ValidationError"


class InternalFailure(AWSError):
    TYPE = "InternalFailure"
    STATUS = 500


class ClientException(AWSError):
    TYPE = "ClientException"
    STATUS = 400
