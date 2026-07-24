from _typeshed import Incomplete
from datetime import date, datetime
from enum import Enum
from typing import Final, Literal

from braintree.address import Address
from braintree.credit_card_verification import CreditCardVerification
from braintree.error_result import ErrorResult
from braintree.resource import Resource
from braintree.resource_collection import ResourceCollection
from braintree.subscription import Subscription
from braintree.successful_result import SuccessfulResult

class CreditCard(Resource):
    class CardType:
        AmEx: Final = "American Express"
        CarteBlanche: Final = "Carte Blanche"
        ChinaUnionPay: Final = "China UnionPay"
        DinersClubInternational: Final = "Diners Club"
        Discover: Final = "Discover"
        Electron: Final = "Electron"
        Elo: Final = "Elo"
        Hiper: Final = "Hiper"
        Hipercard: Final = "Hipercard"
        JCB: Final = "JCB"
        Laser: Final = "Laser"
        UK_Maestro: Final = "UK Maestro"
        Maestro: Final = "Maestro"
        MasterCard: Final = "MasterCard"
        Solo: Final = "Solo"
        Switch: Final = "Switch"
        Visa: Final = "Visa"
        Unknown: Final = "Unknown"

    class CustomerLocation:
        International: Final = "international"
        US: Final = "us"

    class CardTypeIndicator:
        Yes: Final = "Yes"
        No: Final = "No"
        Unknown: Final = "Unknown"

    class DebitNetwork(Enum):
        Accel = "ACCEL"
        Maestro = "MAESTRO"
        Nyce = "NYCE"
        Pulse = "PULSE"
        Star = "STAR"
        Star_Access = "STAR_ACCESS"

    Commercial: type[CardTypeIndicator]
    DurbinRegulated: type[CardTypeIndicator]
    Debit: type[CardTypeIndicator]
    Healthcare: type[CardTypeIndicator]
    CountryOfIssuance: type[CardTypeIndicator]
    IssuingBank: type[CardTypeIndicator]
    Payroll: type[CardTypeIndicator]
    Prepaid: type[CardTypeIndicator]
    ProductId: type[CardTypeIndicator]
    PrepaidReloadable: type[CardTypeIndicator]
    Business: type[CardTypeIndicator]
    Consumer: type[CardTypeIndicator]
    Corporate: type[CardTypeIndicator]
    Purchase: type[CardTypeIndicator]
    @staticmethod
    def create(params: dict[str, Incomplete] | None = None) -> SuccessfulResult | ErrorResult | None: ...
    @staticmethod
    def update(credit_card_token: str, params: dict[str, Incomplete] | None = None) -> SuccessfulResult | ErrorResult | None: ...
    @staticmethod
    def delete(credit_card_token: str) -> SuccessfulResult: ...
    @staticmethod
    def expired() -> ResourceCollection: ...
    @staticmethod
    def expiring_between(start_date: date | datetime, end_date: date | datetime) -> ResourceCollection: ...
    @staticmethod
    def find(credit_card_token: str) -> CreditCard: ...
    @staticmethod
    def from_nonce(nonce: str) -> CreditCard: ...
    @staticmethod
    def create_signature() -> list[str | dict[str, list[str]] | dict[str, list[str | dict[str, list[str]]]]]: ...
    @staticmethod
    def update_signature() -> list[str | dict[str, list[str]] | dict[str, list[str | dict[str, list[str]]]]]: ...
    @staticmethod
    def signature(
        type: Literal["create", "update", "update_via_customer"],
    ) -> list[str | dict[str, list[str]] | dict[str, list[str | dict[str, list[str]]]]]: ...
    is_expired = expired
    billing_address: Address | None
    subscriptions: list[Subscription]
    verification: CreditCardVerification
    def __init__(self, gateway, attributes) -> None: ...
    @property
    def expiration_date(self) -> str | None: ...
    @property
    def masked_number(self) -> str: ...
