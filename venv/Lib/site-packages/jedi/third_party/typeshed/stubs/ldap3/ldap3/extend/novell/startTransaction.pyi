from _typeshed import Incomplete

from ...extend.operation import ExtendedOperation
from ...protocol.novell import CreateGroupTypeRequestValue, CreateGroupTypeResponseValue

class StartTransaction(ExtendedOperation):
    request_name: str
    response_name: str
    request_value: CreateGroupTypeRequestValue
    asn1_spec: CreateGroupTypeResponseValue
    def config(self) -> None: ...
    def __init__(self, connection, controls=None) -> None: ...
    def populate_result(self) -> None: ...
    response_value: Incomplete
    def set_response(self) -> None: ...
