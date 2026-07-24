from ...extend.operation import ExtendedOperation
from ...protocol.novell import Identity

class GetBindDn(ExtendedOperation):
    request_name: str
    response_name: str
    response_attribute: str
    asn1_spec: Identity
    def config(self) -> None: ...
    def populate_result(self) -> None: ...
