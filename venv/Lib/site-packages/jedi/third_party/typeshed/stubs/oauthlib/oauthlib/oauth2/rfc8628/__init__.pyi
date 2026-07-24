from logging import Logger

from .errors import (
    AuthorizationPendingError as AuthorizationPendingError,
    ExpiredTokenError as ExpiredTokenError,
    SlowDownError as SlowDownError,
)

log: Logger
