from functools import partial

from moto.utilities.utils import PARTITION_NAMES, get_partition


def make_arn(
    name: str, account_id: str, region_name: str, _id: str, scope: str, resource: str
) -> str:
    if scope == "REGIONAL":
        scope = "regional"
    elif scope == "CLOUDFRONT":
        scope = "global"

    if region_name in PARTITION_NAMES:
        region_name = "global"
    partition = (
        region_name if region_name in PARTITION_NAMES else get_partition(region_name)
    )

    return f"arn:{partition}:wafv2:{region_name}:{account_id}:{scope}/{resource}/{name}/{_id}"


make_arn_for_wacl = partial(make_arn, resource="webacl")
make_arn_for_ip_set = partial(make_arn, resource="ipset")
make_arn_for_logging_configuration = partial(make_arn, resource="loggingconfiguration")
make_arn_for_rule_group = partial(make_arn, resource="rulegroup")
