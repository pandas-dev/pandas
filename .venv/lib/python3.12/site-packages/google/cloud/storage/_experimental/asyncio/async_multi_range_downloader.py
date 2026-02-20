import warnings

# Import everything from the new stable module
from google.cloud.storage.asyncio.async_multi_range_downloader import *  # noqa

warnings.warn(
    "google.cloud.storage._experimental.asyncio.async_multi_range_downloader has been moved to google.cloud.storage.asyncio.async_multi_range_downloader. "
    "Please update your imports.",
    DeprecationWarning,
    stacklevel=2,
)
