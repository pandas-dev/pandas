from typing import Optional

from moto.core.exceptions import JsonRESTError


class AliasExistsException(JsonRESTError):
    def __init__(self) -> None:
        super().__init__(
            "AliasExistsException", "An account with the given email already exists."
        )


class ResourceNotFoundError(JsonRESTError):
    def __init__(self, message: Optional[str]):
        super().__init__(error_type="ResourceNotFoundException", message=message or "")


class UserNotFoundError(JsonRESTError):
    def __init__(self, message: str):
        super().__init__(error_type="UserNotFoundException", message=message)


class UsernameExistsException(JsonRESTError):
    def __init__(self, message: str):
        super().__init__(error_type="UsernameExistsException", message=message)


class GroupExistsException(JsonRESTError):
    def __init__(self, message: str):
        super().__init__(error_type="GroupExistsException", message=message)


class NotAuthorizedError(JsonRESTError):
    def __init__(self, message: Optional[str]):
        super().__init__(error_type="NotAuthorizedException", message=message or "")


class UserNotConfirmedException(JsonRESTError):
    def __init__(self, message: str):
        super().__init__(error_type="UserNotConfirmedException", message=message)


class ExpiredCodeException(JsonRESTError):
    def __init__(self, message: str):
        super().__init__(error_type="ExpiredCodeException", message=message)


class InvalidParameterException(JsonRESTError):
    def __init__(self, msg: Optional[str] = None):
        self.code = 400
        super().__init__(
            "InvalidParameterException", msg or "A parameter is specified incorrectly."
        )


class InvalidPasswordException(JsonRESTError):
    def __init__(self) -> None:
        super().__init__(
            error_type="InvalidPasswordException",
            message="The provided password does not confirm to the configured password policy",
        )
