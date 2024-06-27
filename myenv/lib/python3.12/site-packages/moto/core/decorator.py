from typing import TYPE_CHECKING, Callable, Optional, TypeVar, Union, overload

from moto import settings
from moto.core.config import DefaultConfig
from moto.core.models import MockAWS, ProxyModeMockAWS, ServerModeMockAWS

if TYPE_CHECKING:
    from typing_extensions import ParamSpec

    P = ParamSpec("P")

T = TypeVar("T")


@overload
def mock_aws(func: "Callable[P, T]") -> "Callable[P, T]": ...


@overload
def mock_aws(
    func: None = None, config: Optional[DefaultConfig] = None
) -> "MockAWS": ...


def mock_aws(
    func: "Optional[Callable[P, T]]" = None,
    config: Optional[DefaultConfig] = None,
) -> Union["MockAWS", "Callable[P, T]"]:
    clss = (
        ServerModeMockAWS
        if settings.TEST_SERVER_MODE
        else (ProxyModeMockAWS if settings.is_test_proxy_mode() else MockAWS)
    )
    if func is not None:
        return clss().__call__(func=func)
    else:
        return clss(config)
