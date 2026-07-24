from moto.utilities.utils import get_partition


def make_arn_for_lexicon(account_id: str, name: str, region_name: str) -> str:
    return f"arn:{get_partition(region_name)}:polly:{region_name}:{account_id}:lexicon/{name}"
