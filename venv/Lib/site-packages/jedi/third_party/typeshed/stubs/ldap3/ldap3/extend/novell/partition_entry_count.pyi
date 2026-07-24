from ...protocol.rfc4511 import LDAPDN
from ..operation import ExtendedOperation

class PartitionEntryCount(ExtendedOperation):
    request_name: str
    response_name: str
    request_value: LDAPDN
    response_attribute: str
    def config(self) -> None: ...
    def __init__(self, connection, partition_dn, controls=None) -> None: ...
    def populate_result(self) -> None: ...
