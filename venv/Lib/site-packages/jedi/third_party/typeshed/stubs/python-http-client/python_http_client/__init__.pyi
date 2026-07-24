from typing import Final

from .client import Client as Client
from .exceptions import (
    BadRequestsError as BadRequestsError,
    ForbiddenError as ForbiddenError,
    GatewayTimeoutError as GatewayTimeoutError,
    HTTPError as HTTPError,
    InternalServerError as InternalServerError,
    MethodNotAllowedError as MethodNotAllowedError,
    NotFoundError as NotFoundError,
    PayloadTooLargeError as PayloadTooLargeError,
    ServiceUnavailableError as ServiceUnavailableError,
    TooManyRequestsError as TooManyRequestsError,
    UnauthorizedError as UnauthorizedError,
    UnsupportedMediaTypeError as UnsupportedMediaTypeError,
)

dir_path: Final[str]
__version__: Final[str]
