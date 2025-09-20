from moto.utilities.utils import get_partition


def make_arn_for_load_balancer(
    account_id: str, lb_id: str, name: str, region_name: str
) -> str:
    return f"arn:{get_partition(region_name)}:elasticloadbalancing:{region_name}:{account_id}:loadbalancer/app/{name}/{lb_id}"


def make_arn_for_target_group(
    account_id: str, tg_id: str, name: str, region_name: str
) -> str:
    return f"arn:{get_partition(region_name)}:elasticloadbalancing:{region_name}:{account_id}:targetgroup/{name}/{tg_id}"
