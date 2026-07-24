from moto.moto_api._internal import mock_random


def get_random_identity_id(region: str) -> str:
    return f"{region}:{mock_random.uuid4()}"
