from typing import Final

from braintree.attribute_getter import AttributeGetter
from braintree.error_result import ErrorResult
from braintree.resource_collection import ResourceCollection
from braintree.successful_result import SuccessfulResult
from braintree.us_bank_account import UsBankAccount

class UsBankAccountVerification(AttributeGetter):
    class Status:
        Failed: Final = "failed"
        GatewayRejected: Final = "gateway_rejected"
        ProcessorDeclined: Final = "processor_declined"
        Unrecognized: Final = "unrecognized"
        Verified: Final = "verified"
        Pending: Final = "pending"

    class VerificationMethod:
        NetworkCheck: Final = "network_check"
        IndependentCheck: Final = "independent_check"
        InstantVerificationAccountValidation: Final = "instant_verification_account_validation"
        TokenizedCheck: Final = "tokenized_check"
        MicroTransfers: Final = "micro_transfers"

    class VerificationAddOns:
        CustomerVerification: Final = "customer_verification"

    us_bank_account: UsBankAccount | None
    def __init__(self, gateway, attributes) -> None: ...
    @staticmethod
    def confirm_micro_transfer_amounts(verification_id: str, amounts) -> SuccessfulResult | ErrorResult | None: ...
    @staticmethod
    def find(verification_id: str) -> UsBankAccountVerification: ...
    @staticmethod
    def search(*query) -> ResourceCollection: ...
    def __eq__(self, other: object) -> bool: ...
