from moto.core.exceptions import JsonRESTError


class GreengrassClientError(JsonRESTError):
    code = 400


class IdNotFoundException(GreengrassClientError):
    def __init__(self, msg: str):
        self.code = 404
        super().__init__("IdNotFoundException", msg)


class InvalidContainerDefinitionException(GreengrassClientError):
    def __init__(self, msg: str):
        self.code = 400
        super().__init__("InvalidContainerDefinitionException", msg)


class VersionNotFoundException(GreengrassClientError):
    def __init__(self, msg: str):
        self.code = 404
        super().__init__("VersionNotFoundException", msg)


class InvalidInputException(GreengrassClientError):
    def __init__(self, msg: str):
        self.code = 400
        super().__init__("InvalidInputException", msg)


class MissingCoreException(GreengrassClientError):
    def __init__(self, msg: str):
        self.code = 400
        super().__init__("MissingCoreException", msg)


class ResourceNotFoundException(GreengrassClientError):
    def __init__(self, msg: str):
        self.code = 404
        super().__init__("ResourceNotFoundException", msg)
