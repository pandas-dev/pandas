import sys

if sys.version_info >= (3, 11):
    from enum import StrEnum

    class ExchangeType(StrEnum):
        direct = "direct"
        fanout = "fanout"
        headers = "headers"
        topic = "topic"

else:
    from enum import Enum

    class ExchangeType(str, Enum):
        direct = "direct"
        fanout = "fanout"
        headers = "headers"
        topic = "topic"
