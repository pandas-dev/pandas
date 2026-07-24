from _typeshed import Incomplete
from typing import Final

from braintree.error_result import ErrorResult
from braintree.resource import Resource
from braintree.successful_result import SuccessfulResult

class DocumentUpload(Resource):
    class Kind:
        EvidenceDocument: Final = "evidence_document"

    @staticmethod
    def create(params: dict[str, Incomplete] | None = None) -> SuccessfulResult | ErrorResult: ...
    @staticmethod
    def create_signature() -> list[str]: ...
    def __init__(self, gateway, attributes) -> None: ...
