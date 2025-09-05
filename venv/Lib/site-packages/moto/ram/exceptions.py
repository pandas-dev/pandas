from moto.core.exceptions import JsonRESTError


class InvalidParameterException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidParameterException", message)


class MalformedArnException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("MalformedArnException", message)


class OperationNotPermittedException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "OperationNotPermittedException",
            "Unable to enable sharing with AWS Organizations. "
            "Received AccessDeniedException from AWSOrganizations with the following error message: "
            "You don't have permissions to access this resource.",
        )


class UnknownResourceException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("UnknownResourceException", message)
