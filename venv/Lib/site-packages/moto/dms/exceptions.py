from moto.core.exceptions import JsonRESTError


class DmsClientError(JsonRESTError):
    code = 400


class ResourceNotFoundFault(DmsClientError):
    def __init__(self, message: str):
        super().__init__("ResourceNotFoundFault", message)


class InvalidResourceStateFault(DmsClientError):
    def __init__(self, message: str):
        super().__init__("InvalidResourceStateFault", message)


class ResourceAlreadyExistsFault(DmsClientError):
    def __init__(self, message: str):
        super().__init__("ResourceAlreadyExistsFault", message)
