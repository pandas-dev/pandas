from moto.utilities.utils import get_partition


def make_arn_for_dashboard(account_id: str, region_name: str, name: str) -> str:
    return f"arn:{get_partition(region_name)}:cloudwatch::{account_id}:dashboard/{name}"


def make_arn_for_alarm(region: str, account_id: str, alarm_name: str) -> str:
    return f"arn:{get_partition(region)}:cloudwatch:{region}:{account_id}:alarm:{alarm_name}"
