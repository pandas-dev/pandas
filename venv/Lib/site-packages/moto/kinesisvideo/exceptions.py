from moto.core.exceptions import ServiceException


class KinesisvideoClientError(ServiceException):
    pass


class ResourceNotFoundException(KinesisvideoClientError):
    def __init__(self) -> None:
        super().__init__(
            "ResourceNotFoundException",
            "The requested stream is not found or not active.",
        )


class ResourceInUseException(KinesisvideoClientError):
    def __init__(self, message: str):
        super().__init__("ResourceInUseException", message)
