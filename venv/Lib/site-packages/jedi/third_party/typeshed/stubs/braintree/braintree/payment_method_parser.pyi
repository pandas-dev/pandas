from braintree.amex_express_checkout_card import AmexExpressCheckoutCard
from braintree.android_pay_card import AndroidPayCard
from braintree.apple_pay_card import ApplePayCard
from braintree.credit_card import CreditCard
from braintree.europe_bank_account import EuropeBankAccount
from braintree.masterpass_card import MasterpassCard
from braintree.paypal_account import PayPalAccount
from braintree.samsung_pay_card import SamsungPayCard
from braintree.sepa_direct_debit_account import SepaDirectDebitAccount
from braintree.unknown_payment_method import UnknownPaymentMethod
from braintree.us_bank_account import UsBankAccount
from braintree.venmo_account import VenmoAccount
from braintree.visa_checkout_card import VisaCheckoutCard

def parse_payment_method(
    gateway, attributes
) -> (
    PayPalAccount
    | CreditCard
    | EuropeBankAccount
    | ApplePayCard
    | AndroidPayCard
    | AmexExpressCheckoutCard
    | SepaDirectDebitAccount
    | VenmoAccount
    | UsBankAccount
    | VisaCheckoutCard
    | MasterpassCard
    | SamsungPayCard
    | UnknownPaymentMethod
): ...
