from braintree.attribute_getter import AttributeGetter
from braintree.receiver import Receiver
from braintree.sender import Sender

class Transfer(AttributeGetter):
    sender: Sender
    receiver: Receiver
    def __init__(self, attributes) -> None: ...
