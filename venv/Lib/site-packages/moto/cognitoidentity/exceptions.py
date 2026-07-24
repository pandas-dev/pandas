from moto.core.exceptions import JsonRESTError


class ResourceNotFoundError(JsonRESTError):
    def __init__(self, message: str):
        super().__init__(error_type="ResourceNotFoundException", message=message)


class InvalidNameException(JsonRESTError):
    message = "1 validation error detected: Value '{}' at 'identityPoolName' failed to satisfy constraint: Member must satisfy regular expression pattern: [\\w\\s+=,.@-]+"

    def __init__(self, name: str):
        msg = InvalidNameException.message.format(name)
        super().__init__(error_type="ValidationException", message=msg)
