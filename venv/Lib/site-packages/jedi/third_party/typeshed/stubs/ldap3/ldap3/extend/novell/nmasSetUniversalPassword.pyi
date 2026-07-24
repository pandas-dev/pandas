from ...extend.operation import ExtendedOperation
from ...protocol.novell import NmasSetUniversalPasswordRequestValue, NmasSetUniversalPasswordResponseValue

class NmasSetUniversalPassword(ExtendedOperation):
    request_name: str
    response_name: str
    request_value: NmasSetUniversalPasswordRequestValue
    asn1_spec: NmasSetUniversalPasswordResponseValue
    response_attribute: str
    def config(self) -> None: ...
    def __init__(self, connection, user, new_password, controls=None) -> None: ...
    def populate_result(self) -> None: ...
