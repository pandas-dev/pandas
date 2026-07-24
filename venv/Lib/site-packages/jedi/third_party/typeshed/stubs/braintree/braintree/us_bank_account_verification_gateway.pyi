from _typeshed import Incomplete

from braintree.error_result import ErrorResult
from braintree.resource_collection import ResourceCollection
from braintree.successful_result import SuccessfulResult
from braintree.us_bank_account_verification import UsBankAccountVerification

class UsBankAccountVerificationGateway:
    gateway: Incomplete
    config: Incomplete
    def __init__(self, gateway) -> None: ...
    def confirm_micro_transfer_amounts(self, verification_id: str, amounts) -> SuccessfulResult | ErrorResult | None: ...
    def find(self, verification_id: str) -> UsBankAccountVerification: ...
    def search(self, *query) -> ResourceCollection: ...
