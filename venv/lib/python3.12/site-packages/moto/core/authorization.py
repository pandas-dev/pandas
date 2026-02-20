import functools
from collections.abc import Callable, Generator
from contextlib import contextmanager
from typing import (
    TYPE_CHECKING,
    ClassVar,
    Optional,
    TypeVar,
)
from urllib.parse import urlparse

import requests

from moto import settings
from moto.utilities.utils import get_partition

if TYPE_CHECKING:
    from typing_extensions import ParamSpec

    P = ParamSpec("P")

T = TypeVar("T")


class ActionAuthenticatorMixin:
    request_count: ClassVar[int] = 0

    PUBLIC_OPERATIONS = [
        "AWSCognitoIdentityService.GetId",
        "AWSCognitoIdentityService.GetOpenIdToken",
        "AWSCognitoIdentityProviderService.ConfirmSignUp",
        "AWSCognitoIdentityProviderService.GetUser",
        "AWSCognitoIdentityProviderService.ForgotPassword",
        "AWSCognitoIdentityProviderService.InitiateAuth",
        "AWSCognitoIdentityProviderService.SignUp",
    ]

    def _authenticate_and_authorize_action(
        self, iam_request_cls: type, resource: str = "*"
    ) -> None:
        if (
            ActionAuthenticatorMixin.request_count
            >= settings.INITIAL_NO_AUTH_ACTION_COUNT
        ):
            if (
                self.headers.get("X-Amz-Target")  # type: ignore[attr-defined]
                in ActionAuthenticatorMixin.PUBLIC_OPERATIONS
            ):
                return
            parsed_url = urlparse(self.uri)  # type: ignore[attr-defined]
            path = parsed_url.path
            if parsed_url.query:
                path += "?" + parsed_url.query
            iam_request = iam_request_cls(
                account_id=self.current_account,  # type: ignore[attr-defined]
                method=self.method,  # type: ignore[attr-defined]
                path=path,
                data=self.data,  # type: ignore[attr-defined]
                body=self.raw_body,  # type: ignore[attr-defined]
                headers=self.headers,  # type: ignore[attr-defined]
                action=self._get_action(),  # type: ignore[attr-defined]
            )
            iam_request.check_signature()
            iam_request.check_action_permitted(resource)
        else:
            ActionAuthenticatorMixin.request_count += 1

    def _authenticate_and_authorize_normal_action(self, resource: str = "*") -> None:
        from moto.iam.access_control import IAMRequest

        self._authenticate_and_authorize_action(IAMRequest, resource)

    def _authenticate_and_authorize_s3_action(
        self, bucket_name: Optional[str] = None, key_name: Optional[str] = None
    ) -> None:
        arn = f"{bucket_name or '*'}/{key_name}" if key_name else (bucket_name or "*")
        resource = f"arn:{get_partition(self.region)}:s3:::{arn}"  # type: ignore[attr-defined]

        from moto.iam.access_control import S3IAMRequest

        self._authenticate_and_authorize_action(S3IAMRequest, resource)

    @staticmethod
    def set_initial_no_auth_action_count(
        initial_no_auth_action_count: int,
    ) -> "Callable[[Callable[P, T]], Callable[P, T]]":
        _test_server_mode_endpoint = settings.test_server_mode_endpoint()

        def decorator(function: "Callable[P, T]") -> "Callable[P, T]":
            def wrapper(*args: "P.args", **kwargs: "P.kwargs") -> T:
                if settings.TEST_SERVER_MODE:
                    response = requests.post(
                        f"{_test_server_mode_endpoint}/moto-api/reset-auth",
                        data=str(initial_no_auth_action_count).encode("utf-8"),
                    )
                    original_initial_no_auth_action_count = response.json()[
                        "PREVIOUS_INITIAL_NO_AUTH_ACTION_COUNT"
                    ]
                else:
                    original_initial_no_auth_action_count = (
                        settings.INITIAL_NO_AUTH_ACTION_COUNT
                    )
                    original_request_count = ActionAuthenticatorMixin.request_count
                    settings.INITIAL_NO_AUTH_ACTION_COUNT = initial_no_auth_action_count
                    ActionAuthenticatorMixin.request_count = 0
                try:
                    result = function(*args, **kwargs)
                finally:
                    if settings.TEST_SERVER_MODE:
                        requests.post(
                            f"{_test_server_mode_endpoint}/moto-api/reset-auth",
                            data=str(original_initial_no_auth_action_count).encode(
                                "utf-8"
                            ),
                        )
                    else:
                        ActionAuthenticatorMixin.request_count = original_request_count
                        settings.INITIAL_NO_AUTH_ACTION_COUNT = (
                            original_initial_no_auth_action_count
                        )
                return result

            functools.update_wrapper(wrapper, function)
            wrapper.__wrapped__ = function  # type: ignore[attr-defined]
            return wrapper

        return decorator


@contextmanager
def enable_iam_authentication() -> Generator[None, None, None]:
    """Make it so that all fixtures and tests run with a simulation of IAM permissions."""

    old_initial_no_auth_action_count = settings.INITIAL_NO_AUTH_ACTION_COUNT
    old_request_count = ActionAuthenticatorMixin.request_count
    settings.INITIAL_NO_AUTH_ACTION_COUNT = 0
    ActionAuthenticatorMixin.request_count = 0
    try:
        yield
    finally:
        settings.INITIAL_NO_AUTH_ACTION_COUNT = old_initial_no_auth_action_count
        ActionAuthenticatorMixin.request_count = old_request_count


@contextmanager
def disable_iam_authentication() -> Generator[None, None, None]:
    """Inverse of enable()."""

    old_initial_no_auth_action_count = settings.INITIAL_NO_AUTH_ACTION_COUNT
    old_request_count = ActionAuthenticatorMixin.request_count
    settings.INITIAL_NO_AUTH_ACTION_COUNT = float("inf")
    ActionAuthenticatorMixin.request_count = 0
    try:
        yield
    finally:
        settings.INITIAL_NO_AUTH_ACTION_COUNT = old_initial_no_auth_action_count
        ActionAuthenticatorMixin.request_count = old_request_count
