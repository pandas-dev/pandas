from ...extend.operation import ExtendedOperation
from ...protocol.novell import NmasGetUniversalPasswordRequestValue, NmasGetUniversalPasswordResponseValue

class NmasGetUniversalPassword(ExtendedOperation):
    request_name: str
    response_name: str
    request_value: NmasGetUniversalPasswordRequestValue
    asn1_spec: NmasGetUniversalPasswordResponseValue
    response_attribute: str
    def config(self) -> None: ...
    def __init__(self, connection, user, controls=None) -> None: ...
    def populate_result(self) -> None: ...
