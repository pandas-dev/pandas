from enum import Enum

VALID_AUTH_MODE_KEYS = ["Type", "Passwords"]
VALID_ENGINE_TYPES = ["redis", "valkey"]


class AuthenticationTypes(str, Enum):
    NOPASSWORD = "no-password-required"
    PASSWORD = "password"
    IAM = "iam"
