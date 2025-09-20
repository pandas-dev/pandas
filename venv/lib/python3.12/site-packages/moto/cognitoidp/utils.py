import base64
import hashlib
import hmac
import re
import string
from typing import Any, Dict, List, Optional

from moto.moto_api._internal import mock_random as random

FORMATS = {
    "email": r"\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b",
    "phone_number": r"\+\d{,15}",
}


PAGINATION_MODEL = {
    "list_user_pools": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 60,
        "unique_attribute": "arn",
    },
    "list_user_pool_clients": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 60,
        "unique_attribute": "id",
    },
    "list_identity_providers": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 60,
        "unique_attribute": "name",
    },
    "list_users": {
        "input_token": "pagination_token",
        "limit_key": "limit",
        "limit_default": 60,
        "unique_attribute": "id",
    },
    "list_groups": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 60,
        "unique_attribute": "group_name",
    },
    "admin_list_groups_for_user": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 60,
        "unique_attribute": "group_name",
    },
    "list_users_in_group": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 60,
        "unique_attribute": "id",
    },
    "list_resource_servers": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 60,
        "unique_attribute": "identifier",
    },
}


def create_id() -> str:
    size = 26
    chars = list(range(10)) + list(string.ascii_lowercase)
    return "".join(str(random.choice(chars)) for x in range(size))


def check_secret_hash(
    app_client_secret: str,
    app_client_id: str,
    username: str,
    secret_hash: Optional[str],
) -> bool:
    key = bytes(str(app_client_secret).encode("latin-1"))
    msg = bytes(str(username + app_client_id).encode("latin-1"))
    new_digest = hmac.new(key, msg, hashlib.sha256).digest()
    SECRET_HASH = base64.b64encode(new_digest).decode()
    return SECRET_HASH == secret_hash


def validate_username_format(username: str, _format: str = "email") -> bool:
    # if the value of the `_format` param other than `email` or `phone_number`,
    # the default value for the regex will match nothing and the
    # method will return None
    return re.fullmatch(FORMATS.get(_format, r"a^"), username) is not None


def flatten_attrs(attrs: List[Dict[str, Any]]) -> Dict[str, Any]:
    return {attr["Name"]: attr["Value"] for attr in attrs}


def expand_attrs(attrs: Dict[str, Any]) -> List[Dict[str, Any]]:
    return [{"Name": k, "Value": v} for k, v in attrs.items()]


ID_HASH_STRATEGY = "HASH"


def generate_id(strategy: Optional[str], *args: Any) -> str:
    if strategy == ID_HASH_STRATEGY:
        return _generate_id_hash(args)
    else:
        return _generate_id_uuid()


def _generate_id_uuid() -> str:
    return random.uuid4().hex


def _generate_id_hash(args: Any) -> str:
    hasher = hashlib.sha256()

    for arg in args:
        hasher.update(str(arg).encode())

    return hasher.hexdigest()
