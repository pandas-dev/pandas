from moto.core.exceptions import JsonRESTError


class NotFoundException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("NotFoundException", message)
