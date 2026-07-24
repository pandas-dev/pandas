from _typeshed import Incomplete

from pyasn1.type.base import Asn1Type

class ExtendedOperation:
    connection: Incomplete
    decoded_response: Incomplete | None
    result: Incomplete | None
    asn1_spec: Asn1Type | None
    request_name: Incomplete | None
    response_name: Incomplete | None
    request_value: Asn1Type | None
    response_value: Incomplete | None
    response_attribute: Incomplete | None
    controls: Incomplete
    def __init__(self, connection, controls=None) -> None: ...
    def send(self): ...
    def populate_result(self) -> None: ...
    def decode_response(self, response=None) -> None: ...
    def set_response(self) -> None: ...
    def config(self) -> None: ...
