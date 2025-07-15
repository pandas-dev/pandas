from moto.core.exceptions import JsonRESTError


class ResourceNotFoundException(JsonRESTError):
    def __init__(self, message: str):
        super().__init__(error_type="ResourceNotFoundException", message=message)


class ValidationException(JsonRESTError):
    def __init__(self, message: str):
        super().__init__(error_type="ResourceNotFoundException", message=message)
