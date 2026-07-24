from _typeshed import Incomplete
from datetime import datetime
from decimal import Decimal
from typing import Final

from braintree.add_on import AddOn
from braintree.address import Address
from braintree.amex_express_checkout_card import AmexExpressCheckoutCard
from braintree.android_pay_card import AndroidPayCard
from braintree.apple_pay_card import ApplePayCard
from braintree.authorization_adjustment import AuthorizationAdjustment
from braintree.credit_card import CreditCard
from braintree.customer import Customer
from braintree.descriptor import Descriptor
from braintree.disbursement_detail import DisbursementDetail
from braintree.discount import Discount
from braintree.dispute import Dispute
from braintree.europe_bank_account import EuropeBankAccount
from braintree.facilitated_details import FacilitatedDetails
from braintree.facilitator_details import FacilitatorDetails
from braintree.local_payment import LocalPayment
from braintree.masterpass_card import MasterpassCard
from braintree.meta_checkout_card import MetaCheckoutCard
from braintree.meta_checkout_token import MetaCheckoutToken
from braintree.package_details import PackageDetails
from braintree.payment_facilitator import PaymentFacilitator
from braintree.paypal_account import PayPalAccount
from braintree.paypal_here import PayPalHere
from braintree.resource import Resource
from braintree.resource_collection import ResourceCollection
from braintree.risk_data import RiskData
from braintree.samsung_pay_card import SamsungPayCard
from braintree.sepa_direct_debit_account import SepaDirectDebitAccount
from braintree.status_event import StatusEvent
from braintree.subscription_details import SubscriptionDetails
from braintree.three_d_secure_info import ThreeDSecureInfo
from braintree.transfer import Transfer
from braintree.us_bank_account import UsBankAccount
from braintree.venmo_account import VenmoAccount
from braintree.visa_checkout_card import VisaCheckoutCard

class Transaction(Resource):
    class CreatedUsing:
        FullInformation: Final = "full_information"
        Token: Final = "token"

    class GatewayRejectionReason:
        ApplicationIncomplete: Final = "application_incomplete"
        Avs: Final = "avs"
        AvsAndCvv: Final = "avs_and_cvv"
        Cvv: Final = "cvv"
        Duplicate: Final = "duplicate"
        ExcessiveRetry: Final = "excessive_retry"
        Fraud: Final = "fraud"
        RiskThreshold: Final = "risk_threshold"
        ThreeDSecure: Final = "three_d_secure"
        TokenIssuance: Final = "token_issuance"

    class ReasonCode:
        ANY_REASON_CODE: Final = "any_reason_code"

    class Source:
        Api: Final = "api"
        ControlPanel: Final = "control_panel"
        Recurring: Final = "recurring"

    class Status:
        AuthorizationExpired: Final = "authorization_expired"
        Authorized: Final = "authorized"
        Authorizing: Final = "authorizing"
        Failed: Final = "failed"
        GatewayRejected: Final = "gateway_rejected"
        ProcessorDeclined: Final = "processor_declined"
        Settled: Final = "settled"
        SettlementConfirmed: Final = "settlement_confirmed"
        SettlementDeclined: Final = "settlement_declined"
        SettlementFailed: Final = "settlement_failed"
        SettlementPending: Final = "settlement_pending"
        Settling: Final = "settling"
        SubmittedForSettlement: Final = "submitted_for_settlement"
        Voided: Final = "voided"

    class Type:
        Credit: Final = "credit"
        Sale: Final = "sale"

    class IndustryType:
        Lodging: Final = "lodging"
        TravelAndCruise: Final = "travel_cruise"
        TravelAndFlight: Final = "travel_flight"

    class AdditionalCharge:
        Restaurant: Final = "restaurant"
        GiftShop: Final = "gift_shop"
        MiniBar: Final = "mini_bar"
        Telephone: Final = "telephone"
        Laundry: Final = "laundry"
        Other: Final = "other"

    @staticmethod
    def adjust_authorization(transaction_id, amount): ...
    @staticmethod
    def clone_transaction(transaction_id, params): ...
    @staticmethod
    def credit(params=None): ...
    @staticmethod
    def find(transaction_id: str) -> Transaction: ...
    @staticmethod
    def refund(transaction_id, amount_or_options=None): ...
    @staticmethod
    def sale(params=None): ...
    @staticmethod
    def search(*query) -> ResourceCollection: ...
    @staticmethod
    def submit_for_settlement(transaction_id, amount=None, params=None): ...
    @staticmethod
    def update_details(transaction_id, params=None): ...
    @staticmethod
    def void(transaction_id): ...
    @staticmethod
    def create(params): ...
    @staticmethod
    def clone_signature(): ...
    @staticmethod
    def create_signature(): ...
    @staticmethod
    def submit_for_settlement_signature(): ...
    @staticmethod
    def submit_for_partial_settlement_signature(): ...
    @staticmethod
    def package_tracking_signature(): ...
    @staticmethod
    def package_tracking(transaction_id, params=None): ...
    @staticmethod
    def update_details_signature(): ...
    @staticmethod
    def refund_signature(): ...
    @staticmethod
    def submit_for_partial_settlement(transaction_id, amount, params=None): ...
    amount: Decimal
    tax_amount: Decimal | None
    discount_amount: Decimal | None
    shipping_amount: Decimal | None
    billing_details: Address
    credit_card_details: CreditCard
    packages: list[PackageDetails]
    paypal_details: PayPalAccount
    paypal_here_details: PayPalHere
    local_payment_details: LocalPayment
    sepa_direct_debit_account_details: SepaDirectDebitAccount
    europe_bank_account_details: EuropeBankAccount
    us_bank_account: UsBankAccount
    apple_pay_details: ApplePayCard
    android_pay_card_details: AndroidPayCard
    amex_express_checkout_card_details: AmexExpressCheckoutCard
    venmo_account_details: VenmoAccount
    visa_checkout_card_details: VisaCheckoutCard
    masterpass_card_details: MasterpassCard
    samsung_pay_card_details: SamsungPayCard
    meta_checkout_card_details: MetaCheckoutCard
    meta_checkout_token_details: MetaCheckoutToken
    sca_exemption_requested: Incomplete
    customer_details: Customer
    shipping_details: Address
    add_ons: list[AddOn]
    discounts: list[Discount]
    status_history: list[StatusEvent]
    subscription_details: SubscriptionDetails
    descriptor: Descriptor
    disbursement_details: DisbursementDetail
    disputes: list[Dispute]
    authorization_adjustments: list[AuthorizationAdjustment]
    payment_instrument_type: Incomplete
    risk_data: RiskData | None
    three_d_secure_info: ThreeDSecureInfo | None
    facilitated_details: FacilitatedDetails
    facilitator_details: FacilitatorDetails
    network_transaction_id: Incomplete
    payment_facilitator: PaymentFacilitator
    transfer: Transfer
    partially_authorized: bool
    subscription_id: str
    created_at: datetime
    def __init__(self, gateway, attributes) -> None: ...
    @property
    def vault_billing_address(self): ...
    @property
    def vault_credit_card(self): ...
    @property
    def vault_customer(self): ...
    @property
    def is_disbursed(self): ...
    @property
    def line_items(self): ...
