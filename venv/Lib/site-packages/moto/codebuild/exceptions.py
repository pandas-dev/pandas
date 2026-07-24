from moto.core.exceptions import JsonRESTError

""" will need exceptions for each api endpoint hit """


class InvalidInputException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidInputException", message)


class ResourceNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)


class ResourceAlreadyExistsException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ResourceAlreadyExistsException", message)
