"""Exceptions raised by the account service."""

from moto.core.exceptions import JsonRESTError


class UnknownContactType(JsonRESTError):
    def __init__(self, user_arn: str, action: str):
        message = f"User: {user_arn} is not authorized to perform: account:{action} (You specified an invalid Alternate Contact type.)"
        super().__init__(error_type="AccessDeniedException", message=message)


class UnspecifiedContactType(JsonRESTError):
    def __init__(self) -> None:
        super().__init__(
            error_type="ResourceNotFoundException",
            message="No contact of the inputted alternate contact type found.",
        )
