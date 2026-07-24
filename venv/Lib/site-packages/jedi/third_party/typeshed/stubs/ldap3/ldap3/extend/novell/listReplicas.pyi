from ...extend.operation import ExtendedOperation
from ...protocol.novell import ReplicaList
from ...protocol.rfc4511 import LDAPDN

class ListReplicas(ExtendedOperation):
    request_name: str
    response_name: str
    request_value: LDAPDN
    asn1_spec: ReplicaList
    response_attribute: str
    def config(self) -> None: ...
    def __init__(self, connection, server_dn, controls=None) -> None: ...
    def populate_result(self) -> None: ...
