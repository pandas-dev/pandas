from typing import List

from moto.core.exceptions import JsonRESTError


class DomainLimitExceededException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "DomainLimitExceeded",
            "The number of registered domains has exceeded the allowed threshold for this account. If you want to "
            "register more domains please request a higher quota",
        )


class DuplicateRequestException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "DuplicateRequest", "The request is already in progress for the domain."
        )


class InvalidInputException(JsonRESTError):
    code = 400

    def __init__(self, error_msgs: List[str]):
        error_msgs_str = "\n\t".join(error_msgs)
        super().__init__(
            "InvalidInput", f"The requested item is not acceptable.\n\t{error_msgs_str}"
        )


class UnsupportedTLDException(JsonRESTError):
    code = 400

    def __init__(self, tld: str):
        super().__init__(
            "UnsupportedTLD",
            f"Amazon Route53 does not support the top-level domain (TLD) `.{tld}`.",
        )
