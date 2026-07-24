from _typeshed import Incomplete
from typing import Final

from braintree.account_updater_daily_report import AccountUpdaterDailyReport
from braintree.connected_merchant_paypal_status_changed import ConnectedMerchantPayPalStatusChanged
from braintree.connected_merchant_status_transitioned import ConnectedMerchantStatusTransitioned
from braintree.disbursement import Disbursement
from braintree.dispute import Dispute
from braintree.granted_payment_instrument_update import GrantedPaymentInstrumentUpdate
from braintree.local_payment_completed import LocalPaymentCompleted
from braintree.local_payment_expired import LocalPaymentExpired
from braintree.local_payment_funded import LocalPaymentFunded
from braintree.local_payment_reversed import LocalPaymentReversed
from braintree.merchant_account import MerchantAccount
from braintree.oauth_access_revocation import OAuthAccessRevocation
from braintree.partner_merchant import PartnerMerchant
from braintree.payment_method_customer_data_updated_metadata import PaymentMethodCustomerDataUpdatedMetadata
from braintree.resource import Resource
from braintree.revoked_payment_method_metadata import RevokedPaymentMethodMetadata
from braintree.subscription import Subscription
from braintree.transaction import Transaction
from braintree.transaction_review import TransactionReview
from braintree.validation_error_collection import ValidationErrorCollection

class WebhookNotification(Resource):
    class Kind:
        AccountUpdaterDailyReport: Final = "account_updater_daily_report"
        Check: Final = "check"
        ConnectedMerchantPayPalStatusChanged: Final = "connected_merchant_paypal_status_changed"
        ConnectedMerchantStatusTransitioned: Final = "connected_merchant_status_transitioned"
        Disbursement: Final = "disbursement"
        DisbursementException: Final = "disbursement_exception"
        DisputeAccepted: Final = "dispute_accepted"
        DisputeAutoAccepted: Final = "dispute_auto_accepted"
        DisputeDisputed: Final = "dispute_disputed"
        DisputeExpired: Final = "dispute_expired"
        DisputeLost: Final = "dispute_lost"
        DisputeOpened: Final = "dispute_opened"
        DisputeUnderReview: Final = "dispute_under_review"
        DisputeWon: Final = "dispute_won"
        GrantedPaymentMethodRevoked: Final = "granted_payment_method_revoked"
        GrantorUpdatedGrantedPaymentMethod: Final = "grantor_updated_granted_payment_method"
        LocalPaymentCompleted: Final = "local_payment_completed"
        LocalPaymentExpired: Final = "local_payment_expired"
        LocalPaymentFunded: Final = "local_payment_funded"
        LocalPaymentReversed: Final = "local_payment_reversed"
        OAuthAccessRevoked: Final = "oauth_access_revoked"
        PartnerMerchantConnected: Final = "partner_merchant_connected"
        PartnerMerchantDeclined: Final = "partner_merchant_declined"
        PartnerMerchantDisconnected: Final = "partner_merchant_disconnected"
        PaymentMethodCustomerDataUpdated: Final = "payment_method_customer_data_updated"
        PaymentMethodRevokedByCustomer: Final = "payment_method_revoked_by_customer"
        RecipientUpdatedGrantedPaymentMethod: Final = "recipient_updated_granted_payment_method"
        RefundFailed: Final = "refund_failed"
        SubscriptionBillingSkipped: Final = "subscription_billing_skipped"
        SubscriptionCanceled: Final = "subscription_canceled"
        SubscriptionChargedSuccessfully: Final = "subscription_charged_successfully"
        SubscriptionChargedUnsuccessfully: Final = "subscription_charged_unsuccessfully"
        SubscriptionExpired: Final = "subscription_expired"
        SubscriptionTrialEnded: Final = "subscription_trial_ended"
        SubscriptionWentActive: Final = "subscription_went_active"
        SubscriptionWentPastDue: Final = "subscription_went_past_due"
        TransactionDisbursed: Final = "transaction_disbursed"
        TransactionRetried: Final = "transaction_retried"
        TransactionReviewed: Final = "transaction_reviewed"
        TransactionSettled: Final = "transaction_settled"
        TransactionSettlementDeclined: Final = "transaction_settlement_declined"

    @staticmethod
    def parse(signature: str, payload: str) -> WebhookNotification: ...
    @staticmethod
    def verify(challenge: str) -> str: ...
    source_merchant_id: Incomplete
    subscription: Subscription
    merchant_account: MerchantAccount
    transaction: Transaction
    transaction_review: TransactionReview
    connected_merchant_status_transitioned: ConnectedMerchantStatusTransitioned
    connected_merchant_paypal_status_changed: ConnectedMerchantPayPalStatusChanged
    partner_merchant: PartnerMerchant
    oauth_access_revocation: OAuthAccessRevocation
    disbursement: Disbursement
    dispute: Dispute
    account_updater_daily_report: AccountUpdaterDailyReport
    granted_payment_instrument_update: GrantedPaymentInstrumentUpdate
    revoked_payment_method_metadata: RevokedPaymentMethodMetadata
    local_payment_completed: LocalPaymentCompleted
    local_payment_expired: LocalPaymentExpired
    local_payment_funded: LocalPaymentFunded
    local_payment_reversed: LocalPaymentReversed
    payment_method_customer_data_updated_metadata: PaymentMethodCustomerDataUpdatedMetadata
    errors: ValidationErrorCollection
    message: Incomplete
    def __init__(self, gateway, attributes) -> None: ...
