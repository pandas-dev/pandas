from _typeshed import Incomplete

from braintree.error_result import ErrorResult
from braintree.successful_result import SuccessfulResult

class DocumentUploadGateway:
    gateway: Incomplete
    config: Incomplete
    def __init__(self, gateway) -> None: ...
    def create(self, params: dict[str, Incomplete] | None = None) -> SuccessfulResult | ErrorResult: ...
