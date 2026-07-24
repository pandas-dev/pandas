from _typeshed import Incomplete
from decimal import Decimal
from typing import Final

from braintree.attribute_getter import AttributeGetter
from braintree.error_result import ErrorResult
from braintree.resource_collection import ResourceCollection
from braintree.risk_data import RiskData
from braintree.successful_result import SuccessfulResult
from braintree.three_d_secure_info import ThreeDSecureInfo

class CreditCardVerification(AttributeGetter):
    class Status:
        Failed: Final = "failed"
        GatewayRejected: Final = "gateway_rejected"
        ProcessorDeclined: Final = "processor_declined"
        Verified: Final = "verified"

    amount: Decimal | None
    currency_iso_code: Incomplete
    processor_response_code: Incomplete
    processor_response_text: Incomplete
    network_response_code: Incomplete
    network_response_text: Incomplete
    risk_data: RiskData | None
    three_d_secure_info: ThreeDSecureInfo | None
    def __init__(self, gateway, attributes) -> None: ...
    @staticmethod
    def find(verification_id: str) -> CreditCardVerification: ...
    @staticmethod
    def search(*query) -> ResourceCollection: ...
    @staticmethod
    def create(params) -> SuccessfulResult | ErrorResult | None: ...
    @staticmethod
    def create_signature() -> list[dict[str, list[str | dict[str, list[str]]]] | dict[str, list[str]] | str]: ...
    def __eq__(self, other: object) -> bool: ...
