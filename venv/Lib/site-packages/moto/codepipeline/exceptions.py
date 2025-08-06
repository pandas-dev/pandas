from moto.core.exceptions import JsonRESTError


class InvalidStructureException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidStructureException", message)


class PipelineNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("PipelineNotFoundException", message)


class ResourceNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)


class InvalidTagsException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidTagsException", message)


class TooManyTagsException(JsonRESTError):
    code = 400

    def __init__(self, arn: str):
        super().__init__(
            "TooManyTagsException", f"Tag limit exceeded for resource [{arn}]."
        )
