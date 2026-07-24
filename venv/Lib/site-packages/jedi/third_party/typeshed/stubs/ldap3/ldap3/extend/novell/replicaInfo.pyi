from ...protocol.novell import ReplicaInfoRequestValue
from ..operation import ExtendedOperation

class ReplicaInfo(ExtendedOperation):
    request_name: str
    response_name: str
    request_value: ReplicaInfoRequestValue
    response_attribute: str
    def config(self) -> None: ...
    def __init__(self, connection, server_dn, partition_dn, controls=None) -> None: ...
    def populate_result(self) -> None: ...
