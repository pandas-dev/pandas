from ...extend.operation import ExtendedOperation
from ...protocol.rfc3062 import PasswdModifyRequestValue, PasswdModifyResponseValue

class ModifyPassword(ExtendedOperation):
    request_name: str
    request_value: PasswdModifyRequestValue
    asn1_spec: PasswdModifyResponseValue
    response_attribute: str
    def config(self) -> None: ...
    def __init__(
        self, connection, user=None, old_password=None, new_password=None, hash_algorithm=None, salt=None, controls=None
    ) -> None: ...
    def populate_result(self) -> None: ...
