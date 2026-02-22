from moto.core.exceptions import JsonRESTError


class DataBrewClientError(JsonRESTError):
    code = 400


class AlreadyExistsException(DataBrewClientError):
    def __init__(self, typ: str):
        super().__init__("AlreadyExistsException", f"{typ} already exists.")


class ConflictException(DataBrewClientError):
    code = 409

    def __init__(self, message: str):
        super().__init__("ConflictException", message)


class ValidationException(DataBrewClientError):
    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class RulesetAlreadyExistsException(AlreadyExistsException):
    def __init__(self) -> None:
        super().__init__("Ruleset")


class EntityNotFoundException(DataBrewClientError):
    def __init__(self, msg: str):
        super().__init__("EntityNotFoundException", msg)


class ResourceNotFoundException(DataBrewClientError):
    code = 404

    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)


class RulesetNotFoundException(EntityNotFoundException):
    def __init__(self, recipe_name: str):
        super().__init__(f"Ruleset {recipe_name} not found.")


class ServiceQuotaExceededException(JsonRESTError):
    code = 402

    def __init__(self) -> None:
        super().__init__(
            "ServiceQuotaExceededException", "A service quota is exceeded."
        )
