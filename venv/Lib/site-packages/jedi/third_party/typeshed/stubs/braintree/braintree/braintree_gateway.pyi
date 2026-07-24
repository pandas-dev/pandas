from braintree.add_on_gateway import AddOnGateway
from braintree.address_gateway import AddressGateway
from braintree.apple_pay_gateway import ApplePayGateway
from braintree.bank_account_instant_verification_gateway import BankAccountInstantVerificationGateway
from braintree.client_token_gateway import ClientTokenGateway
from braintree.configuration import Configuration
from braintree.credit_card_gateway import CreditCardGateway
from braintree.credit_card_verification_gateway import CreditCardVerificationGateway
from braintree.customer_gateway import CustomerGateway
from braintree.discount_gateway import DiscountGateway
from braintree.dispute_gateway import DisputeGateway
from braintree.document_upload_gateway import DocumentUploadGateway
from braintree.exchange_rate_quote_gateway import ExchangeRateQuoteGateway
from braintree.merchant_account_gateway import MerchantAccountGateway
from braintree.merchant_gateway import MerchantGateway
from braintree.oauth_gateway import OAuthGateway
from braintree.payment_method_gateway import PaymentMethodGateway
from braintree.payment_method_nonce_gateway import PaymentMethodNonceGateway
from braintree.paypal_account_gateway import PayPalAccountGateway
from braintree.paypal_payment_resource_gateway import PayPalPaymentResourceGateway
from braintree.plan_gateway import PlanGateway
from braintree.sepa_direct_debit_account_gateway import SepaDirectDebitAccountGateway
from braintree.settlement_batch_summary_gateway import SettlementBatchSummaryGateway
from braintree.subscription_gateway import SubscriptionGateway
from braintree.testing_gateway import TestingGateway
from braintree.transaction_gateway import TransactionGateway
from braintree.transaction_line_item_gateway import TransactionLineItemGateway
from braintree.us_bank_account_gateway import UsBankAccountGateway
from braintree.us_bank_account_verification_gateway import UsBankAccountVerificationGateway
from braintree.util.graphql_client import GraphQLClient
from braintree.webhook_notification_gateway import WebhookNotificationGateway
from braintree.webhook_testing_gateway import WebhookTestingGateway

class BraintreeGateway:
    config: Configuration
    add_on: AddOnGateway
    address: AddressGateway
    apple_pay: ApplePayGateway
    bank_account_instant_verification: BankAccountInstantVerificationGateway
    client_token: ClientTokenGateway
    credit_card: CreditCardGateway
    customer: CustomerGateway
    discount: DiscountGateway
    dispute: DisputeGateway
    document_upload: DocumentUploadGateway
    exchange_rate_quote: ExchangeRateQuoteGateway
    graphql_client: GraphQLClient
    merchant: MerchantGateway
    merchant_account: MerchantAccountGateway
    oauth: OAuthGateway
    payment_method: PaymentMethodGateway
    payment_method_nonce: PaymentMethodNonceGateway
    paypal_account: PayPalAccountGateway
    paypal_payment_resource: PayPalPaymentResourceGateway
    plan: PlanGateway
    sepa_direct_debit_account: SepaDirectDebitAccountGateway
    settlement_batch_summary: SettlementBatchSummaryGateway
    subscription: SubscriptionGateway
    testing: TestingGateway
    transaction: TransactionGateway
    transaction_line_item: TransactionLineItemGateway
    us_bank_account: UsBankAccountGateway
    us_bank_account_verification: UsBankAccountVerificationGateway
    verification: CreditCardVerificationGateway
    webhook_notification: WebhookNotificationGateway
    webhook_testing: WebhookTestingGateway
    def __init__(self, config=None, **kwargs) -> None: ...
    def close(self) -> None: ...
