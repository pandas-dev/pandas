from moto.core.exceptions import ServiceException


class OrganizationsClientError(ServiceException):
    pass


class AccountAlreadyRegisteredException(OrganizationsClientError):
    def __init__(self) -> None:
        super().__init__(
            "AccountAlreadyRegisteredException",
            "The provided account is already a delegated administrator for your organization.",
        )


class AccountAlreadyClosedException(OrganizationsClientError):
    def __init__(self) -> None:
        super().__init__(
            "AccountAlreadyClosedException",
            "The provided account is already closed.",
        )


class AccountNotRegisteredException(OrganizationsClientError):
    def __init__(self) -> None:
        super().__init__(
            "AccountNotRegisteredException",
            "The provided account is not a registered delegated administrator for your organization.",
        )


class AccountNotFoundException(OrganizationsClientError):
    def __init__(self) -> None:
        super().__init__(
            "AccountNotFoundException", "You specified an account that doesn't exist."
        )


class AlreadyInOrganizationException(OrganizationsClientError):
    def __init__(self) -> None:
        super().__init__(
            "AlreadyInOrganizationException",
            "The provided account is already a member of an organization.",
        )


class AWSOrganizationsNotInUseException(OrganizationsClientError):
    def __init__(self) -> None:
        super().__init__(
            "AWSOrganizationsNotInUseException",
            "Your account is not a member of an organization.",
        )


class ConstraintViolationException(OrganizationsClientError):
    def __init__(self, message: str):
        super().__init__("ConstraintViolationException", message)


class InvalidInputException(OrganizationsClientError):
    def __init__(self, message: str):
        super().__init__("InvalidInputException", message)


class DuplicateOrganizationalUnitException(OrganizationsClientError):
    def __init__(self) -> None:
        super().__init__(
            "DuplicateOrganizationalUnitException",
            "An OU with the same name already exists.",
        )


class DuplicatePolicyException(OrganizationsClientError):
    def __init__(self) -> None:
        super().__init__(
            "DuplicatePolicyException", "A policy with the same name already exists."
        )


class OrganizationNotEmptyException(OrganizationsClientError):
    def __init__(self) -> None:
        super().__init__(
            "OrganizationNotEmptyException",
            "To delete an organization you must first remove all member accounts (except the master).",
        )


class PolicyTypeAlreadyEnabledException(OrganizationsClientError):
    def __init__(self) -> None:
        super().__init__(
            "PolicyTypeAlreadyEnabledException",
            "The specified policy type is already enabled.",
        )


class PolicyTypeNotEnabledException(OrganizationsClientError):
    def __init__(self) -> None:
        super().__init__(
            "PolicyTypeNotEnabledException",
            "This operation can be performed only for enabled policy types.",
        )


class RootNotFoundException(OrganizationsClientError):
    def __init__(self) -> None:
        super().__init__(
            "RootNotFoundException", "You specified a root that doesn't exist."
        )


class TargetNotFoundException(OrganizationsClientError):
    def __init__(self) -> None:
        super().__init__(
            "TargetNotFoundException", "You specified a target that doesn't exist."
        )


class PolicyNotFoundException(OrganizationsClientError):
    def __init__(self, message: str) -> None:
        super().__init__("PolicyNotFoundException", message)
