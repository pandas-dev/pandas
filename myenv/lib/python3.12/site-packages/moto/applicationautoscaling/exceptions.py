from moto.core.exceptions import JsonRESTError


class AWSValidationException(JsonRESTError):
    def __init__(self, message: str) -> None:
        super().__init__("ValidationException", message)
