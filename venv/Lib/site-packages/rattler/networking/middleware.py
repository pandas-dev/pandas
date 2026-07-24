from __future__ import annotations

from typing import Callable

from rattler.rattler import (
    PyAddHeadersMiddleware,
    PyAuthenticationMiddleware,
    PyGCSMiddleware,
    PyMirrorMiddleware,
    PyOciMiddleware,
    PyRetryMiddleware,
    PyS3Middleware,
    PyS3Config,
)


class MirrorMiddleware:
    def __init__(self, mirrors: dict[str, list[str]]) -> None:
        """
        Create a new MirrorMiddleware instance.
        The mirrors argument should be a dictionary where the keys are the
        original mirror URLs and the values are lists of mirror URLs to
        replace the original mirror with.

        Examples
        --------
        ```python
        >>> from rattler.networking import Client
        >>> middleware = MirrorMiddleware({"https://conda.anaconda.org/conda-forge": ["https://repo.prefix.dev/conda-forge"]})
        >>> middleware
        MirrorMiddleware()
        >>> Client([middleware])
        Client()
        >>>
        ```
        """
        self._middleware = PyMirrorMiddleware(mirrors)

    def __repr__(self) -> str:
        """
        Returns a representation of the Middleware

        Examples
        --------
        ```python
        >>> middleware = MirrorMiddleware({"https://conda.anaconda.org/conda-forge": ["https://repo.prefix.dev/conda-forge"]})
        >>> middleware
        MirrorMiddleware()
        >>>
        ```
        """
        return f"{type(self).__name__}()"


class AuthenticationMiddleware:
    """
    Middleware to handle authentication from keychain
    """

    def __init__(self) -> None:
        self._middleware = PyAuthenticationMiddleware()

    def __repr__(self) -> str:
        """
        Returns a representation of the Middleware

        Examples
        --------
        ```python
        >>> from rattler.networking import Client
        >>> middleware = AuthenticationMiddleware()
        >>> middleware
        AuthenticationMiddleware()
        >>> Client([middleware])
        Client()
        >>>
        ```
        """
        return f"{type(self).__name__}()"


class RetryMiddleware:
    """
    Middleware that retries transient HTTP errors with exponential back-off.

    Examples
    --------
    ```python
    >>> from rattler.networking import Client, RetryMiddleware
    >>> RetryMiddleware()
    RetryMiddleware(max_retries=3)
    >>> RetryMiddleware(5)
    RetryMiddleware(max_retries=5)
    >>> Client([RetryMiddleware()])
    Client()
    >>>
    ```
    """

    def __init__(self, max_retries: int = 3) -> None:
        self._max_retries = max_retries
        self._middleware = PyRetryMiddleware(max_retries)

    def __repr__(self) -> str:
        return f"{type(self).__name__}(max_retries={self._max_retries})"


class OciMiddleware:
    """
    Middleware to handle `oci://` URLs
    """

    def __init__(self) -> None:
        self._middleware = PyOciMiddleware()

    def __repr__(self) -> str:
        """
        Returns a representation of the Middleware

        Examples
        --------
        ```python
        >>> from rattler.networking import Client
        >>> middleware = OciMiddleware()
        >>> middleware
        OciMiddleware()
        >>> Client([middleware])
        Client()
        >>>
        ```
        """
        return f"{type(self).__name__}()"


class GCSMiddleware:
    """
    Middleware to work with gcs:// URLs

    Examples
    --------
    ```python
    >>> from rattler.networking import Client
    >>> middleware = GCSMiddleware()
    >>> middleware
    GCSMiddleware()
    >>> Client([middleware])
    Client()
    >>>
    ```
    """

    def __init__(self) -> None:
        self._middleware = PyGCSMiddleware()

    def __repr__(self) -> str:
        return f"{type(self).__name__}()"


class S3Config:
    """
    Middleware to work with s3:// URLs

    Examples
    --------
    ```python
    >>> from rattler.networking import S3Middleware
    >>> config = S3Config("http://localhost:9000", "eu-central-1", True)
    >>> config
    S3Config(http://localhost:9000, eu-central-1, True)
    >>> middleware = S3Middleware({"my-bucket": config})
    >>> middleware
    S3Middleware()
    >>> S3Config()
    S3Config(aws sdk)
    >>>
    ```
    """

    def __init__(
        self, endpoint_url: str | None = None, region: str | None = None, force_path_style: bool | None = None
    ) -> None:
        self._config = PyS3Config(endpoint_url, region, force_path_style)
        if (endpoint_url is None) != (region is None) or (endpoint_url is None) != (force_path_style is None):
            raise ValueError("Invalid arguments for S3Config")
        self._endpoint_url = endpoint_url
        self._region = region
        self._force_path_style = force_path_style

    def __repr__(self) -> str:
        inner = (
            f"{self._endpoint_url}, {self._region}, {self._force_path_style}"
            if self._endpoint_url is not None
            else "aws sdk"
        )
        return f"{type(self).__name__}({inner})"


class S3Middleware:
    """
    Middleware to work with s3:// URLs

    Examples
    --------
    ```python
    >>> from rattler.networking import Client
    >>> middleware = S3Middleware()
    >>> middleware
    S3Middleware()
    >>> Client([middleware])
    Client()
    >>>
    ```
    """

    def __init__(self, config: dict[str, S3Config] | None = None) -> None:
        if config is None:
            config = dict()
        self._middleware = PyS3Middleware({k: v._config for k, v in config.items()})

    def __repr__(self) -> str:
        return f"{type(self).__name__}()"


class AddHeadersMiddleware:
    """
    Middleware that adds headers to requests based on a callback function.

    The callback receives the host and path of the request URL and should return
    a dictionary of headers to add, or None to add no headers.

    Examples
    --------
    ```python
    >>> from rattler.networking import Client, AddHeadersMiddleware
    >>> def my_callback(host: str, path: str) -> dict[str, str] | None:
    ...     if host == "private.example.com":
    ...         return {"Authorization": "Bearer my-token"}
    ...     return None
    >>> middleware = AddHeadersMiddleware(my_callback)
    >>> middleware
    AddHeadersMiddleware()
    >>> Client([middleware])
    Client()
    >>>
    ```
    """

    def __init__(self, callback: Callable[[str, str], dict[str, str] | None]) -> None:
        """
        Create a new AddHeadersMiddleware instance.

        Args:
            callback: A callable that takes (host, path) and returns a dictionary
                of headers to add to the request, or None to add no headers.
                The host is the hostname of the request URL (e.g., "conda.anaconda.org").
                The path is the path component of the URL (e.g., "/conda-forge/linux-64/repodata.json").
        """
        self._middleware = PyAddHeadersMiddleware(callback)

    def __repr__(self) -> str:
        return f"{type(self).__name__}()"
