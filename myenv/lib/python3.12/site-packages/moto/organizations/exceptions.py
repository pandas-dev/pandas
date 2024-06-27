from moto.core.exceptions import JsonRESTError


class AccountAlreadyRegisteredException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "AccountAlreadyRegisteredException",
            "The provided account is already a delegated administrator for your organization.",
        )


class AccountNotRegisteredException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "AccountNotRegisteredException",
            "The provided account is not a registered delegated administrator for your organization.",
        )


class AccountNotFoundException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "AccountNotFoundException", "You specified an account that doesn't exist."
        )


class AWSOrganizationsNotInUseException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "AWSOrganizationsNotInUseException",
            "Your account is not a member of an organization.",
        )


class ConstraintViolationException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ConstraintViolationException", message)


class InvalidInputException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidInputException", message)


class DuplicateOrganizationalUnitException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "DuplicateOrganizationalUnitException",
            "An OU with the same name already exists.",
        )


class DuplicatePolicyException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "DuplicatePolicyException", "A policy with the same name already exists."
        )


class PolicyTypeAlreadyEnabledException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "PolicyTypeAlreadyEnabledException",
            "The specified policy type is already enabled.",
        )


class PolicyTypeNotEnabledException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "PolicyTypeNotEnabledException",
            "This operation can be performed only for enabled policy types.",
        )


class RootNotFoundException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "RootNotFoundException", "You specified a root that doesn't exist."
        )


class TargetNotFoundException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "TargetNotFoundException", "You specified a target that doesn't exist."
        )


class PolicyNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, message: str) -> None:
        super().__init__("PolicyNotFoundException", message)
