from _typeshed import Incomplete

class BaseEndpoint:
    request_validator: Incomplete
    token_generator: Incomplete
    def __init__(self, request_validator, token_generator=None) -> None: ...
