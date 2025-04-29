import base64
import os
import string

from moto.moto_api._internal import mock_random as random

ACCOUNT_SPECIFIC_ACCESS_KEY_PREFIX = "8NWMTLYQ"
ACCOUNT_SPECIFIC_ASSUMED_ROLE_ID_PREFIX = "3X42LBCD"
SESSION_TOKEN_PREFIX = "FQoGZXIvYXdzEBYaD"
DEFAULT_STS_SESSION_DURATION = 3600


def random_session_token() -> str:
    return (
        SESSION_TOKEN_PREFIX
        + base64.b64encode(os.urandom(266))[len(SESSION_TOKEN_PREFIX) :].decode()
    )


def random_assumed_role_id() -> str:
    return (
        ACCOUNT_SPECIFIC_ASSUMED_ROLE_ID_PREFIX + _random_uppercase_or_digit_sequence(9)
    )


def _random_uppercase_or_digit_sequence(length: int) -> str:
    return "".join(
        str(random.choice(string.ascii_uppercase + string.digits))
        for _ in range(length)
    )
