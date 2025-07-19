from moto.core.exceptions import JsonRESTError


class ConflictException(JsonRESTError):
    def __init__(self, message: str):
        super().__init__("ConflictException", message)


class BadRequestException(JsonRESTError):
    def __init__(self, message: str):
        super().__init__("BadRequestException", message)
