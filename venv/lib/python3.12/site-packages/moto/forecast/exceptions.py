from moto.core.exceptions import AWSError


class InvalidInputException(AWSError):
    TYPE = "InvalidInputException"


class ResourceAlreadyExistsException(AWSError):
    TYPE = "ResourceAlreadyExistsException"


class ResourceNotFoundException(AWSError):
    TYPE = "ResourceNotFoundException"


class ResourceInUseException(AWSError):
    TYPE = "ResourceInUseException"


class LimitExceededException(AWSError):
    TYPE = "LimitExceededException"


class ValidationException(AWSError):
    TYPE = "ValidationException"
