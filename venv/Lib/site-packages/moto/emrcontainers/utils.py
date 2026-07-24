import string

from moto.moto_api._internal import mock_random as random


def random_id(size: int = 13) -> str:
    chars = list(range(10)) + list(string.ascii_lowercase)
    return "".join(str(random.choice(chars)) for x in range(size))


def random_cluster_id() -> str:
    return random_id(size=25)


def random_job_id() -> str:
    return random_id(size=19)
