import string
from email.utils import parseaddr
from typing import Optional

from moto.moto_api._internal import mock_random as random


def random_hex(length: int) -> str:
    return "".join(random.choice(string.ascii_lowercase) for x in range(length))


def get_random_message_id() -> str:
    return f"{random_hex(16)}-{random_hex(8)}-{random_hex(4)}-{random_hex(4)}-{random_hex(4)}-{random_hex(12)}-{random_hex(6)}"


def is_valid_address(addr: str) -> Optional[str]:
    _, address = parseaddr(addr)
    address_parts = address.split("@")
    if len(address_parts) != 2 or not address_parts[1]:
        return "Missing domain"
    return None
