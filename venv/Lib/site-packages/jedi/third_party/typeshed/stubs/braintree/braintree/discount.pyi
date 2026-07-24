from braintree.modification import Modification

class Discount(Modification):
    @staticmethod
    def all() -> list[Discount]: ...
