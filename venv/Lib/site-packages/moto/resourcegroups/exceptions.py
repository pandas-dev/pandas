from moto.core.exceptions import JsonRESTError


class BadRequestException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__(error_type="BadRequestException", message=message)


class NotFoundException(JsonRESTError):
    code = 404

    def __init__(
        self, message: str = "One or more of the specified resources don't exist."
    ):
        super().__init__(error_type="NotFoundException", message=message)
