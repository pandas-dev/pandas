"""Exceptions raised by the qldb service."""

from moto.core.exceptions import JsonRESTError


class QLDBError(JsonRESTError):
    code = 400


class LedgerNotFoundException(QLDBError):
    def __init__(self, name: str) -> None:
        super().__init__(
            "LedgerNotFound",
            f"{name} does not match any existing ledger.",
        )


class LedgerNameAlreadyTakenException(QLDBError):
    def __init__(self, name: str) -> None:
        super().__init__(
            "LedgerNameAlreadyTaken",
            f"There is already a ledger with the name '{name}'",
        )


class ResourceNotFoundException(QLDBError):
    def __init__(self, resource_arn: str) -> None:
        super().__init__(
            "ResourceNotFound",
            f"{resource_arn} does not match the resource arn of any existing ledger",
        )
