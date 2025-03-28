import inspect
import re

from boto3 import Session

from moto.eks.exceptions import InvalidParameterException


def get_partition(region: str) -> str:
    valid_matches = [
        # (region prefix, aws partition)
        ("cn-", "aws-cn"),
        ("us-gov-", "aws-us-gov"),
        ("us-gov-iso-", "aws-iso"),
        ("us-gov-iso-b-", "aws-iso-b"),
    ]

    for prefix, partition in valid_matches:
        if region.startswith(prefix):
            return partition
    return "aws"


def method_name(use_parent: bool = False) -> str:
    """
    Returns the name of the method which called it from the stack in PascalCase.
    If `use_parent` is True, returns the parent of the method which called it instead.
    For example:  False/default will return the name of the method calling it.
    In a helper method, use True to return the name of the method which called the helper.
    """
    return (
        # stack()[0] is this method, stack()[1] is the method which called this one, etc
        inspect.stack()[int(use_parent) + 1][0]
        .f_code.co_name.replace("_", " ")
        .title()
        .replace(" ", "")
    )


def validate_role_arn(arn: str) -> None:
    valid_role_arn_format = re.compile(
        "arn:(?P<partition>.+):iam::(?P<account_id>[0-9]{12}):role/.+"
    )
    match = valid_role_arn_format.match(arn)
    valid_partition = match.group("partition") in Session().get_available_partitions()  # type: ignore

    if not all({arn, match, valid_partition}):
        raise InvalidParameterException("Invalid Role Arn: '" + arn + "'")  # type: ignore
