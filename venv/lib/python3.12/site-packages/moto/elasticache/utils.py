from enum import Enum

PAGINATION_MODEL = {
    "describe_cache_clusters": {
        "input_token": "marker",
        "limit_key": "max_records",
        "limit_default": 100,
        "unique_attribute": "cache_cluster_id",
    },
    "describe_cache_subnet_groups": {
        "input_token": "marker",
        "limit_key": "max_records",
        "limit_default": 100,
        "unique_attribute": "cache_subnet_group_name",
    },
    "describe_replication_groups": {
        "input_token": "marker",
        "limit_key": "max_records",
        "limit_default": 100,
        "unique_attribute": "replication_group_id",
    },
}
VALID_AUTH_MODE_KEYS = ["Type", "Passwords"]
VALID_ENGINE_TYPES = ["redis", "valkey"]


class AuthenticationTypes(str, Enum):
    NOPASSWORD = "no-password-required"
    PASSWORD = "password"
    IAM = "iam"
