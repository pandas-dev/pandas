import asyncio
import json
import logging
import random

import aiohttp.client_exceptions
import google.auth.exceptions
import requests.exceptions
from decorator import decorator
from google.api_core import exceptions as api_exceptions
from google.api_core.retry import AsyncRetry

logger = logging.getLogger("gcsfs")
DEFAULT_RETRY_CONFIG = {
    "timeout": 60.0,
    "initial": 1.0,
    "maximum": 60.0,
    "multiplier": 2.0,
}


class HttpError(Exception):
    """Holds the message and code from cloud errors."""

    def __init__(self, error_response=None):
        # Save error_response for potential pickle.
        self._error_response = error_response
        if error_response:
            self.code = error_response.get("code", None)
            self.message = error_response.get("message", "")
            if self.code:
                if isinstance(self.message, bytes):
                    self.message += (", %s" % self.code).encode()
                else:
                    self.message += ", %s" % self.code
        else:
            self.message = ""
            self.code = None
        # Call the base class constructor with the parameters it needs
        super().__init__(self.message)

    def __reduce__(self):
        """This makes the Exception pickleable."""

        # This is basically deconstructing the HttpError when pickled.
        return HttpError, (self._error_response,)


class ChecksumError(Exception):
    """Raised when the md5 hash of the content does not match the header."""

    pass


class NonRetryableError(Exception):
    """Raised when the underlying error can not be retried, or continued further."""

    pass


RETRIABLE_EXCEPTIONS = (
    requests.exceptions.ChunkedEncodingError,
    requests.exceptions.ConnectionError,
    requests.exceptions.ReadTimeout,
    requests.exceptions.Timeout,
    requests.exceptions.ProxyError,
    requests.exceptions.SSLError,
    requests.exceptions.ContentDecodingError,
    google.auth.exceptions.RefreshError,
    aiohttp.client_exceptions.ClientError,
    ChecksumError,
)


errs = list(range(500, 505)) + [
    # Request Timeout
    408,
    # Too Many Requests
    429,
]
errs = set(errs + [str(e) for e in errs])


def is_retriable(exception):
    """Returns True if this exception is retriable."""
    if isinstance(exception, NonRetryableError):
        return False

    if isinstance(exception, HttpError):
        # Add 401 to retriable errors when it's an auth expiration issue
        if exception.code == 401 and "Invalid Credentials" in str(exception.message):
            return True
        return exception.code in errs

    return isinstance(exception, RETRIABLE_EXCEPTIONS)


def validate_response(status, content, path, args=None):
    """
    Check the requests object r, raise error if it's not ok.

    Parameters
    ----------
    r: requests response object
    path: associated URL path, for error messages
    """
    if status >= 400 and status != 499:
        # 499 is special "upload was cancelled" status
        if args:
            from .core import quote

            path = path.format(*[quote(p) for p in args])
        if status == 404:
            raise FileNotFoundError(path)

        error = None
        msg = ""
        if content:
            if hasattr(content, "decode"):
                content = content.decode()
            try:
                error = json.loads(content)["error"]
                # Sometimes the error message is a string.
                if isinstance(error, str):
                    msg = error
                else:
                    msg = error["message"]
            except json.decoder.JSONDecodeError:
                msg = content

        if status == 403:
            raise OSError(f"Forbidden: {path}\n{msg}")
        elif status == 412:
            raise FileExistsError(path)
        elif status == 502:
            raise requests.exceptions.ProxyError()
        elif "invalid" in str(msg):
            raise ValueError(f"Bad Request: {path}\n{msg}")
        elif error and not isinstance(error, str):
            raise HttpError(error)
        elif status:
            raise HttpError({"code": status, "message": msg})  # text-like
        else:
            raise RuntimeError(msg)


@decorator
async def retry_request(func, retries=6, *args, **kwargs):
    for retry in range(retries):
        try:
            if retry > 0:
                await asyncio.sleep(min(random.random() + 2 ** (retry - 1), 32))
            return await func(*args, **kwargs)
        except (
            HttpError,
            requests.exceptions.RequestException,
            google.auth.exceptions.GoogleAuthError,
            ChecksumError,
            aiohttp.client_exceptions.ClientError,
        ) as e:
            if (
                isinstance(e, HttpError)
                and e.code == 400
                and "requester pays" in e.message
            ):
                msg = (
                    "Bucket is requester pays. "
                    "Set `requester_pays=True` when creating the GCSFileSystem."
                )
                raise ValueError(msg) from e
            # Special test for 404 to avoid retrying the request
            if (
                isinstance(e, aiohttp.client_exceptions.ClientResponseError)
                and e.status == 404
            ):
                logger.debug("Request returned 404, no retries.")
                raise e
            if isinstance(e, HttpError) and e.code == 404:
                logger.debug("Request returned 404, no retries.")
                raise e
            if retry == retries - 1:
                logger.exception(f"{func.__name__} out of retries on exception: {e}")
                raise e
            if is_retriable(e):
                logger.debug(f"{func.__name__} retrying after exception: {e}")
                continue
            logger.exception(f"{func.__name__} non-retriable exception: {e}")
            raise e


def _is_transient_exception(exception):
    is_transient = isinstance(
        exception,
        (
            api_exceptions.DeadlineExceeded,
            api_exceptions.ServiceUnavailable,
            api_exceptions.InternalServerError,
            api_exceptions.TooManyRequests,
            api_exceptions.ResourceExhausted,
            api_exceptions.Unknown,
        ),
    )
    if (
        not is_transient
        and isinstance(exception, api_exceptions.Unauthenticated)
        and "Invalid Credentials" in str(exception)
    ):
        is_transient = True
    return is_transient


def get_storage_control_retry_config(base_config=None, **kwargs) -> AsyncRetry:
    """
    Returns an AsyncRetry object configured for Storage Control API calls.

    Priority: kwargs (timeout, etc.) > base_config > package defaults.

    Args:
        base_config: A dict containing base settings.
        **kwargs: Direct call-site overrides (e.g., timeout=10).
    """
    retry_kwargs = DEFAULT_RETRY_CONFIG.copy()
    valid_keys = DEFAULT_RETRY_CONFIG.keys()
    if base_config:
        retry_kwargs.update(
            {k: v for k, v in base_config.items() if k in valid_keys and v is not None}
        )

    overrides = {k: v for k, v in kwargs.items() if k in valid_keys and v is not None}
    retry_kwargs.update(overrides)

    return AsyncRetry(predicate=_is_transient_exception, **retry_kwargs)
