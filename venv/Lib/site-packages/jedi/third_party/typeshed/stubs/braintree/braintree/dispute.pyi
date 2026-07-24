from _typeshed import Incomplete
from decimal import Decimal
from typing import Final

from braintree.attribute_getter import AttributeGetter
from braintree.dispute_details import DisputeEvidence, DisputePayPalMessage, DisputeStatusHistory
from braintree.error_result import ErrorResult
from braintree.successful_result import SuccessfulResult
from braintree.transaction_details import TransactionDetails

class Dispute(AttributeGetter):
    class Status:
        Accepted: Final = "accepted"
        AutoAccepted: Final = "auto_accepted"
        Disputed: Final = "disputed"
        Expired: Final = "expired"
        Lost: Final = "lost"
        Open: Final = "open"
        UnderReview: Final = "under_review"
        Won: Final = "won"

    class Reason:
        CancelledRecurringTransaction: Final = "cancelled_recurring_transaction"
        CreditNotProcessed: Final = "credit_not_processed"
        Duplicate: Final = "duplicate"
        Fraud: Final = "fraud"
        General: Final = "general"
        InvalidAccount: Final = "invalid_account"
        NotRecognized: Final = "not_recognized"
        ProductNotReceived: Final = "product_not_received"
        ProductUnsatisfactory: Final = "product_unsatisfactory"
        Retrieval: Final = "retrieval"
        TransactionAmountDiffers: Final = "transaction_amount_differs"

    class Kind:
        Chargeback: Final = "chargeback"
        PreArbitration: Final = "pre_arbitration"
        Retrieval: Final = "retrieval"

    class ChargebackProtectionLevel:
        Effortless: Final = "effortless"
        Standard: Final = "standard"
        NotProtected: Final = "not_protected"

    class PreDisputeProgram:
        NONE: Final = "none"
        VisaRdr: Final = "visa_rdr"

    class ProtectionLevel:
        EffortlessCBP: Final = "Effortless Chargeback Protection tool"
        StandardCBP: Final = "Chargeback Protection tool"
        NoProtection: Final = "No Protection"

    @staticmethod
    def accept(id: str) -> SuccessfulResult | ErrorResult: ...
    @staticmethod
    def add_file_evidence(dispute_id: str, document_upload_id) -> SuccessfulResult | ErrorResult | None: ...
    @staticmethod
    def add_text_evidence(id: str, content_or_request) -> SuccessfulResult | ErrorResult | None: ...
    @staticmethod
    def finalize(id: str) -> SuccessfulResult | ErrorResult: ...
    @staticmethod
    def find(id: str) -> Dispute: ...
    @staticmethod
    def remove_evidence(id: str, evidence_id: str) -> SuccessfulResult | ErrorResult: ...
    @staticmethod
    def search(*query) -> SuccessfulResult: ...
    amount: Decimal | None
    amount_disputed: Decimal | None
    amount_won: Decimal | None
    protection_level: Incomplete
    transaction_details: TransactionDetails
    transaction = transaction_details
    evidence: list[DisputeEvidence] | None
    paypal_messages: list[DisputePayPalMessage] | None
    status_history: list[DisputeStatusHistory] | None
    processor_comments: Incomplete
    forwarded_comments: processor_comments
    def __init__(self, attributes) -> None: ...
