from moto.core.exceptions import JsonRESTError


class BadRequestException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__(error_type="BadRequestException", message=message)
