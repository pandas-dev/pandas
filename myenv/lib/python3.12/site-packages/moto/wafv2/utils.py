from moto.utilities.utils import PARTITION_NAMES, get_partition


def make_arn_for_wacl(
    name: str, account_id: str, region_name: str, wacl_id: str, scope: str
) -> str:
    """https://docs.aws.amazon.com/waf/latest/developerguide/how-aws-waf-works.html - explains --scope (cloudfront vs regional)"""
    return make_arn(name, account_id, region_name, wacl_id, scope, "webacl")


def make_arn_for_ip_set(
    name: str, account_id: str, region_name: str, _id: str, scope: str
) -> str:
    return make_arn(name, account_id, region_name, _id, scope, "ipset")


def make_arn_for_logging_configuration(
    name: str, account_id: str, region_name: str, _id: str, scope: str
) -> str:
    return make_arn(name, account_id, region_name, _id, scope, "logingconfiguration")


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
