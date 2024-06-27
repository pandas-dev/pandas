from moto.core.exceptions import JsonRESTError


class InvalidRequestException(JsonRESTError):
    def __init__(self, message: str):
        super().__init__("InvalidRequestException", message)


class ValidationException(JsonRESTError):
    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class ResourceNotFoundException(JsonRESTError):
    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)
