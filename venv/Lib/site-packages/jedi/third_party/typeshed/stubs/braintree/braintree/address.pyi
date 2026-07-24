from _typeshed import Incomplete
from typing import Final

from braintree.error_result import ErrorResult
from braintree.resource import Resource
from braintree.successful_result import SuccessfulResult

class Address(Resource):
    class ShippingMethod:
        SameDay: Final = "same_day"
        NextDay: Final = "next_day"
        Priority: Final = "priority"
        Ground: Final = "ground"
        Electronic: Final = "electronic"
        ShipToStore: Final = "ship_to_store"
        PickupInStore: Final = "pickup_in_store"

    @staticmethod
    def create(params: dict[str, Incomplete] | None = None) -> SuccessfulResult | ErrorResult | None: ...
    @staticmethod
    def delete(customer_id: str, address_id: str) -> SuccessfulResult: ...
    @staticmethod
    def find(customer_id: str, address_id: str) -> Address: ...
    @staticmethod
    def update(
        customer_id: str, address_id: str, params: dict[str, Incomplete] | None = None
    ) -> SuccessfulResult | ErrorResult | None: ...
    @staticmethod
    def create_signature() -> list[str | dict[str, list[str]]]: ...
    @staticmethod
    def update_signature() -> list[str | dict[str, list[str]]]: ...
