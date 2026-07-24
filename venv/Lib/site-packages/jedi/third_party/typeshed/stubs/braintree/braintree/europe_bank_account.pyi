from typing import Final

from braintree.resource import Resource

class EuropeBankAccount(Resource):
    class MandateType:
        Business: Final = "business"
        Consumer: Final = "consumer"

    @staticmethod
    def signature() -> list[str]: ...
