from typing import Final

from braintree.resource import Resource

class MerchantAccount(Resource):
    class Status:
        Active: Final = "active"
        Pending: Final = "pending"
        Suspended: Final = "suspended"

    class FundingDestination:
        Bank: Final = "bank"
        Email: Final = "email"
        MobilePhone: Final = "mobile_phone"

    FundingDestinations: type[FundingDestination]
    master_merchant_account: MerchantAccount
    def __init__(self, gateway, attributes) -> None: ...
    @staticmethod
    def create(params=None): ...
    @staticmethod
    def update(id, attributes): ...
    @staticmethod
    def find(id): ...
