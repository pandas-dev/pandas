from _typeshed import Incomplete
from datetime import date
from decimal import Decimal
from typing import Final

from braintree.add_on import AddOn
from braintree.descriptor import Descriptor
from braintree.discount import Discount
from braintree.error_result import ErrorResult
from braintree.resource import Resource
from braintree.resource_collection import ResourceCollection
from braintree.subscription_status_event import SubscriptionStatusEvent
from braintree.successful_result import SuccessfulResult
from braintree.transaction import Transaction

class Subscription(Resource):
    class TrialDurationUnit:
        Day: Final = "day"
        Month: Final = "month"

    class Source:
        Api: Final = "api"
        ControlPanel: Final = "control_panel"
        Recurring: Final = "recurring"

    class Status:
        Active: Final = "Active"
        Canceled: Final = "Canceled"
        Expired: Final = "Expired"
        PastDue: Final = "Past Due"
        Pending: Final = "Pending"

    @staticmethod
    def create(params: dict[str, Incomplete] | None = None) -> SuccessfulResult | ErrorResult | None: ...
    @staticmethod
    def create_signature(): ...
    @staticmethod
    def find(subscription_id: str) -> Subscription: ...
    @staticmethod
    def retry_charge(subscription_id, amount=None, submit_for_settlement: bool = False): ...
    @staticmethod
    def update(subscription_id: str, params: dict[str, Incomplete] | None = None) -> SuccessfulResult | ErrorResult | None: ...
    @staticmethod
    def cancel(subscription_id: str) -> SuccessfulResult | ErrorResult | None: ...
    @staticmethod
    def search(*query) -> ResourceCollection: ...
    @staticmethod
    def update_signature(): ...
    price: Decimal
    balance: Decimal
    next_billing_date: date
    next_billing_period_amount: Decimal
    add_ons: list[AddOn]
    descriptor: Descriptor
    description: Incomplete
    discounts: list[Discount]
    status_history: list[SubscriptionStatusEvent]
    transactions: list[Transaction]
    def __init__(self, gateway, attributes) -> None: ...
