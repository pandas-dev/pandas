import base64
import string

from moto.moto_api._internal import mock_random as random

AWS_ROLE_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ234567"
ACCOUNT_OFFSET = 549755813888  # int.from_bytes(base64.b32decode(b"QAAAAAAA"), byteorder="big"), start value


def _random_uppercase_or_digit_sequence(length: int) -> str:
    return "".join(str(random.choice(AWS_ROLE_ALPHABET)) for _ in range(length))


def generate_access_key_id_from_account_id(
    account_id: str, prefix: str, total_length: int = 20
) -> str:
    """
    Generates a key id (e.g. access key id) for the given account id and prefix

    :param account_id: Account id this key id should belong to
    :param prefix: Prefix, e.g. ASIA for temp credentials or AROA for roles
    :param total_length: Total length of the access key (e.g. 20 for temp access keys, 21 for role ids)
    :return: Generated id
    """
    account_id_nr = int(account_id)
    id_with_offset = account_id_nr // 2 + ACCOUNT_OFFSET
    account_bytes = int.to_bytes(id_with_offset, byteorder="big", length=5)
    account_part = base64.b32encode(account_bytes).decode("utf-8")
    middle_char = (
        random.choice(AWS_ROLE_ALPHABET[16:])
        if account_id_nr % 2
        else random.choice(AWS_ROLE_ALPHABET[:16])
    )
    semi_fixed_part = prefix + account_part + middle_char
    return semi_fixed_part + _random_uppercase_or_digit_sequence(
        total_length - len(semi_fixed_part)
    )


def random_alphanumeric(length: int) -> str:
    options = string.ascii_letters + string.digits + "+" + "/"
    return "".join(random.choices(options, k=length))


def random_resource_id(size: int = 20) -> str:
    return "".join(random.choices(string.ascii_lowercase + string.digits, k=size))


def random_role_id(account_id: str) -> str:
    return generate_access_key_id_from_account_id(
        account_id=account_id, prefix="AROA", total_length=21
    )


def random_access_key() -> str:
    return "".join(random.choices(string.ascii_uppercase + string.digits, k=16))


def random_policy_id() -> str:
    return "A" + "".join(random.choices(string.ascii_uppercase + string.digits, k=20))
