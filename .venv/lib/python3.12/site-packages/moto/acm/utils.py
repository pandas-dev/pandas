from moto.moto_api._internal import mock_random
from moto.utilities.utils import get_partition


def make_arn_for_certificate(account_id: str, region_name: str) -> str:
    # Example
    # arn:aws:acm:eu-west-2:764371465172:certificate/c4b738b8-56fe-4b3a-b841-1c047654780b
    return f"arn:{get_partition(region_name)}:acm:{region_name}:{account_id}:certificate/{mock_random.uuid4()}"
