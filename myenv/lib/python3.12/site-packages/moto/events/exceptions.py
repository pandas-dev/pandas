from typing import Optional

from moto.core.exceptions import JsonRESTError


class IllegalStatusException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("IllegalStatusException", message)


class InvalidEventPatternException(JsonRESTError):
    code = 400

    def __init__(self, reason: Optional[str] = None):
        msg = "Event pattern is not valid. "
        if reason:
            msg += f"Reason: {reason}"

        super().__init__("InvalidEventPatternException", msg)


class ResourceNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)


class ResourceAlreadyExistsException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ResourceAlreadyExistsException", message)


class ValidationException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ValidationException", message)
