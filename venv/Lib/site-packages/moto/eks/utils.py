import re

from boto3 import Session

from moto.eks.exceptions import InvalidParameterException


def validate_role_arn(arn: str) -> None:
    valid_partitions = Session().get_available_partitions()
    valid_role_arn_format = re.compile(
        "arn:(?P<partition>.+):iam::(?P<account_id>[0-9]{12}):role/.+"
    )
    match = valid_role_arn_format.match(arn)
    if not all(
        [
            arn,
            match,
            match.group("partition") in valid_partitions if match else False,
        ]
    ):
        raise InvalidParameterException(message="Invalid Role Arn: '" + arn + "'")
