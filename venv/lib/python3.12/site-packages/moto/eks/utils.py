import re

from boto3 import Session

from moto.eks.exceptions import InvalidParameterException


def validate_role_arn(arn: str) -> None:
    valid_role_arn_format = re.compile(
        "arn:(?P<partition>.+):iam::(?P<account_id>[0-9]{12}):role/.+"
    )
    match = valid_role_arn_format.match(arn)
    valid_partition = match.group("partition") in Session().get_available_partitions()  # type: ignore

    if not all({arn, match, valid_partition}):
        raise InvalidParameterException("Invalid Role Arn: '" + arn + "'")  # type: ignore
