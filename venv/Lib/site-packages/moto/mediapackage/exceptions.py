from moto.core.exceptions import JsonRESTError


class MediaPackageClientError(JsonRESTError):
    code = 400


# AWS service exceptions are caught with the underlying botocore exception, ClientError
class ClientError(MediaPackageClientError):
    def __init__(self, error: str, message: str):
        super().__init__(error, message)
