from _typeshed import Incomplete

from braintree.amex_express_checkout_card import AmexExpressCheckoutCard
from braintree.android_pay_card import AndroidPayCard
from braintree.apple_pay_card import ApplePayCard
from braintree.credit_card import CreditCard
from braintree.error_result import ErrorResult
from braintree.europe_bank_account import EuropeBankAccount
from braintree.masterpass_card import MasterpassCard
from braintree.paypal_account import PayPalAccount
from braintree.samsung_pay_card import SamsungPayCard
from braintree.sepa_direct_debit_account import SepaDirectDebitAccount
from braintree.successful_result import SuccessfulResult
from braintree.unknown_payment_method import UnknownPaymentMethod
from braintree.us_bank_account import UsBankAccount
from braintree.venmo_account import VenmoAccount
from braintree.visa_checkout_card import VisaCheckoutCard

class PaymentMethodGateway:
    gateway: Incomplete
    config: Incomplete
    def __init__(self, gateway) -> None: ...
    def create(self, params: dict[str, Incomplete] | None = None) -> SuccessfulResult | ErrorResult: ...
    def find(
        self, payment_method_token: str
    ) -> (
        AndroidPayCard
        | ApplePayCard
        | EuropeBankAccount
        | CreditCard
        | PayPalAccount
        | UsBankAccount
        | VenmoAccount
        | VisaCheckoutCard
        | AmexExpressCheckoutCard
        | SepaDirectDebitAccount
        | MasterpassCard
        | SamsungPayCard
        | UnknownPaymentMethod
    ): ...
    def update(self, payment_method_token: str, params) -> SuccessfulResult | ErrorResult: ...
    def delete(self, payment_method_token: str, options=None) -> SuccessfulResult: ...
    options: dict[str, Incomplete]
    def grant(self, payment_method_token: str, options=None) -> SuccessfulResult | ErrorResult: ...
    def revoke(self, payment_method_token: str) -> SuccessfulResult | ErrorResult: ...
