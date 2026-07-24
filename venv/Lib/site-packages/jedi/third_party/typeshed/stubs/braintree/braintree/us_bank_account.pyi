from braintree.ach_mandate import AchMandate
from braintree.error_result import ErrorResult
from braintree.resource import Resource
from braintree.successful_result import SuccessfulResult
from braintree.us_bank_account_verification import UsBankAccountVerification

class UsBankAccount(Resource):
    @staticmethod
    def find(token: str) -> UsBankAccount | None: ...
    @staticmethod
    def sale(token: str, transactionRequest) -> SuccessfulResult | ErrorResult | None: ...
    @staticmethod
    def signature() -> list[str]: ...
    ach_mandate: AchMandate | None
    verifications: list[UsBankAccountVerification]
    def __init__(self, gateway, attributes) -> None: ...
