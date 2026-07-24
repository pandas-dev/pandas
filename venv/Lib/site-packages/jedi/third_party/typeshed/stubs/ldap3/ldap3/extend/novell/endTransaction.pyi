from _typeshed import Incomplete

from ...extend.operation import ExtendedOperation
from ...protocol.novell import EndGroupTypeRequestValue, EndGroupTypeResponseValue

class EndTransaction(ExtendedOperation):
    request_name: str
    response_name: str
    request_value: EndGroupTypeRequestValue
    asn1_spec: EndGroupTypeResponseValue
    def config(self) -> None: ...
    def __init__(self, connection, commit: bool = True, controls=None) -> None: ...
    def populate_result(self) -> None: ...
    response_value: Incomplete
    def set_response(self) -> None: ...
