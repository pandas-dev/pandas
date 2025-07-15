from moto.core.exceptions import ServiceException


class InvalidRequestException(ServiceException):
    def __init__(self, message: str):
        super().__init__("InvalidRequestException", message)


class ValidationException(ServiceException):
    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class ResourceNotFoundException(ServiceException):
    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)
