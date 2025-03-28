from moto.core.exceptions import RESTError


class InvalidFormat(RESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__(InvalidFormat.__name__, message)


class InvalidParameterValue(RESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__(InvalidParameterValue.__name__, message)


class InvalidParameterCombination(RESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__(InvalidParameterCombination.__name__, message)


class ResourceNotFound(RESTError):
    code = 404

    def __init__(self) -> None:
        super().__init__(ResourceNotFound.__name__, "Unknown")


class ResourceNotFoundException(RESTError):
    code = 404

    def __init__(self) -> None:
        super().__init__(ResourceNotFoundException.__name__, "Unknown")


class ValidationError(RESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__(ValidationError.__name__, message)
