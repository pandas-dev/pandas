from moto.core.exceptions import ServiceException


class SagemakerClientError(ServiceException):
    pass


class ValidationError(SagemakerClientError):
    code = "ValidationException"

    # We have to keep this initializer because call sites are using message="message"
    def __init__(self, message: str):
        super().__init__(message)


class MissingModel(ValidationError):
    def __init__(self, model: str):
        super().__init__(f'Could not find model "{model}".')


class AWSValidationException(SagemakerClientError):
    code = "ValidationException"


class ResourceInUseException(SagemakerClientError):
    code = "ResourceInUse"

    # We have to keep this initializer because call sites are using message="message"
    def __init__(self, message: str):
        super().__init__(message)


class ConflictException(SagemakerClientError):
    code = "ConflictException"


class ResourceNotFound(SagemakerClientError):
    code = "ResourceNotFound"

    # We have to keep this initializer because call sites are using message="message"
    def __init__(self, message: str):
        super().__init__(message)
