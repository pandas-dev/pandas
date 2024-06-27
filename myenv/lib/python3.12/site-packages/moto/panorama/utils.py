import base64
import hashlib
from datetime import datetime
from typing import Any

from moto.utilities.utils import get_partition


def deep_convert_datetime_to_isoformat(obj: Any) -> Any:
    if isinstance(obj, datetime):
        return obj.isoformat()
    elif isinstance(obj, list):
        return [deep_convert_datetime_to_isoformat(x) for x in obj]
    elif isinstance(obj, dict):
        return {k: deep_convert_datetime_to_isoformat(v) for k, v in obj.items()}
    else:
        return obj


def hash_device_name(name: str) -> str:
    digest = hashlib.md5(name.encode("utf-8")).digest()
    token = base64.b64encode(digest)
    return token.decode("utf-8")


def generate_package_id(name: str) -> str:
    digest = hashlib.md5(name.encode("utf-8")).digest()
    token = base64.b64encode(digest).decode("utf-8")
    # make token only containing alphanumeric characters
    token = "".join(char for char in token if char.isalnum())
    token = token * (64 // len(token) + 1)
    return token[:64].lower()


def arn_formatter(_type: str, _id: str, account_id: str, region_name: str) -> str:
    return f"arn:{get_partition(region_name)}:panorama:{region_name}:{account_id}:{_type}/{_id}"
