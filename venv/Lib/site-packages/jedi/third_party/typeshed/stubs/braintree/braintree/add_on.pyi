from braintree.modification import Modification

class AddOn(Modification):
    @staticmethod
    def all() -> list[AddOn]: ...
