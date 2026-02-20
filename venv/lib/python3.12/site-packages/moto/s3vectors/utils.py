from moto.utilities.utils import get_partition


def create_vector_bucket_arn(account_id: str, region: str, name: str) -> str:
    return f"arn:{get_partition(region)}:s3vectors:{region}:{account_id}:bucket/{name}"
