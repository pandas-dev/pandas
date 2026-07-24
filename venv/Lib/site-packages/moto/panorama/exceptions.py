from moto.core.exceptions import JsonRESTError


class ValidationError(JsonRESTError):
    def __init__(self, message: str):
        super().__init__("ValidationException", message)
