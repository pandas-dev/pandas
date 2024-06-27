from typing import Any

from moto.core.exceptions import RESTError

EXCEPTION_RESPONSE = """<?xml version="1.0"?>
<ErrorResponse xmlns="http://elasticache.amazonaws.com/doc/2015-02-02/">
  <Error>
    <Type>Sender</Type>
    <Code>{{ error_type }}</Code>
    <Message>{{ message }}</Message>
  </Error>
  <{{ request_id_tag }}>30c0dedb-92b1-4e2b-9be4-1188e3ed86ab</{{ request_id_tag }}>
</ErrorResponse>"""


class ElastiCacheException(RESTError):
    code = 400
    extended_templates = {"ecerror": EXCEPTION_RESPONSE}
    env = RESTError.extended_environment(extended_templates)

    def __init__(self, code: str, message: str, **kwargs: Any):
        kwargs.setdefault("template", "ecerror")
        super().__init__(code, message)


class PasswordTooShort(ElastiCacheException):
    code = 404

    def __init__(self) -> None:
        super().__init__(
            "InvalidParameterValue",
            message="Passwords length must be between 16-128 characters.",
        )


class PasswordRequired(ElastiCacheException):
    code = 404

    def __init__(self) -> None:
        super().__init__(
            "InvalidParameterValue",
            message="No password was provided. If you want to create/update the user without password, please use the NoPasswordRequired flag.",
        )


class UserAlreadyExists(ElastiCacheException):
    code = 404

    def __init__(self) -> None:
        super().__init__("UserAlreadyExists", message="User user1 already exists.")


class UserNotFound(ElastiCacheException):
    code = 404

    def __init__(self, user_id: str):
        super().__init__("UserNotFound", message=f"User {user_id} not found.")


class CacheClusterAlreadyExists(ElastiCacheException):
    code = 404

    def __init__(self, cache_cluster_id: str):
        (
            super().__init__(
                "CacheClusterAlreadyExists",
                message=f"Cache cluster {cache_cluster_id} already exists.",
            ),
        )


class CacheClusterNotFound(ElastiCacheException):
    code = 404

    def __init__(self, cache_cluster_id: str):
        super().__init__(
            "CacheClusterNotFound",
            message=f"Cache cluster {cache_cluster_id} not found.",
        )


class InvalidARNFault(ElastiCacheException):
    code = 400

    def __init__(self, arn: str):
        super().__init__(
            "InvalidARNFault",
            message=f"ARN {arn} is invalid.",
        )
