from functools import partial

from moto.utilities.utils import get_partition


def make_arn(
    name: str, account_id: str, region_name: str, _id: str, scope: str, resource: str
) -> str:
    if scope == "REGIONAL":
        scope = "regional"
    elif scope == "CLOUDFRONT":
        scope = "global"
        # cloudfront global scope is managed from us-east-1 region: https://docs.aws.amazon.com/waf/latest/APIReference/API_CreateIPSet.html#WAF-CreateIPSet-request-Scope
        # base_backend stores global region as "aws": https://github.com/getmoto/moto/blob/d00aa025b6c3d37977508b5d5e81ecad4ca15159/moto/core/base_backend.py#L272
        # Therefore needs to be overriden here to correctly form the ARN
        region_name = "us-east-1"

    partition = get_partition(region_name)

    return f"arn:{partition}:wafv2:{region_name}:{account_id}:{scope}/{resource}/{name}/{_id}"


make_arn_for_wacl = partial(make_arn, resource="webacl")
make_arn_for_ip_set = partial(make_arn, resource="ipset")
make_arn_for_logging_configuration = partial(make_arn, resource="loggingconfiguration")
make_arn_for_rule_group = partial(make_arn, resource="rulegroup")


def make_arn_for_regex_pattern_set(
    name: str, account_id: str, region_name: str, _id: str, scope: str
) -> str:
    """
    Generate an ARN for a WAFv2 Regex Pattern Set
    """
    if scope == "CLOUDFRONT":
        region_name = "us-east-1"
    return f"arn:aws:wafv2:{region_name}:{account_id}:regional/regexpatternset/{name}/{_id}"
