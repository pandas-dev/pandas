from moto.core.exceptions import ServiceException


class S3ControlError(ServiceException):
    pass


class AccessPointNotFound(S3ControlError):
    code = "NoSuchAccessPoint"

    def __init__(self, name: str):
        super().__init__("The specified accesspoint does not exist")
        self.access_point_name = name


class AccessPointPolicyNotFound(S3ControlError):
    code = "NoSuchAccessPointPolicy"

    def __init__(self, name: str):
        super().__init__("The specified accesspoint policy does not exist")
        self.access_point_name = name


class MultiRegionAccessPointNotFound(S3ControlError):
    code = "NoSuchMultiRegionAccessPoint"

    def __init__(self, name: str):
        super().__init__("The specified multi-region access point does not exist")
        self.name = name


class MultiRegionAccessPointPolicyNotFound(S3ControlError):
    code = "NoSuchMultiRegionAccessPointPolicy"

    def __init__(self, name: str):
        super().__init__(
            "The specified multi-region access point policy does not exist"
        )
        self.name = name


class MultiRegionAccessPointOperationNotFound(S3ControlError):
    code = "NoSuchAsyncRequest"

    def __init__(self, request_token: str):
        super().__init__("The specified async request does not exist")
        self.request_token_arn = request_token


class NoSuchPublicAccessBlockConfiguration(S3ControlError):
    # Note that this exception is in the different format then the S3 exception with the same name
    # This exception should return a nested response `<ErrorResponse><Error>..`
    # The S3 variant uses a flat `<Error>`-response
    code = "NoSuchPublicAccessBlockConfiguration"

    def __init__(self) -> None:
        super().__init__("The public access block configuration was not found")


class InvalidRequestException(S3ControlError):
    code = "InvalidRequest"

    def __init__(self, message: str):
        super().__init__(message)


class StorageLensConfigurationNotFound(S3ControlError):
    code = "NoSuchConfiguration"

    def __init__(self, config_id: str):
        super().__init__(f"The specified configuration does not exist: {config_id}")
