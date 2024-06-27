import string
from typing import Any, Dict

from moto.moto_api._internal import mock_random as random


def random_id(size: int = 13) -> str:
    chars = list(range(10)) + list(string.ascii_lowercase)
    return "".join(str(random.choice(chars)) for x in range(size))


def random_appplication_id() -> str:
    return random_id(size=16)


def random_job_id() -> str:
    return random_id(size=16)


def default_auto_start_configuration() -> Dict[str, bool]:
    return {"enabled": True}


def default_auto_stop_configuration() -> Dict[str, Any]:
    return {"enabled": True, "idleTimeoutMinutes": 15}
