from ...extend.operation import ExtendedOperation

class WhoAmI(ExtendedOperation):
    request_name: str
    response_attribute: str
    def config(self) -> None: ...
    def populate_result(self) -> None: ...
