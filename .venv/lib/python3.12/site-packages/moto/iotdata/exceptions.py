from moto.core.exceptions import JsonRESTError


class IoTDataPlaneClientError(JsonRESTError):
    code = 400


class ResourceNotFoundException(IoTDataPlaneClientError):
    def __init__(self) -> None:
        self.code = 404
        super().__init__(
            "ResourceNotFoundException", "The specified resource does not exist"
        )


class InvalidRequestException(IoTDataPlaneClientError):
    def __init__(self, message: str):
        self.code = 400
        super().__init__("InvalidRequestException", message)


class ConflictException(IoTDataPlaneClientError):
    def __init__(self, message: str):
        self.code = 409
        super().__init__("ConflictException", message)
