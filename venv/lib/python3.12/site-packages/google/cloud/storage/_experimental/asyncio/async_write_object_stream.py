import warnings

# Import everything from the new stable module
from google.cloud.storage.asyncio.async_write_object_stream import *  # noqa

warnings.warn(
    "google.cloud.storage._experimental.asyncio.async_write_object_stream has been moved to google.cloud.storage.asyncio.async_write_object_stream. "
    "Please update your imports.",
    DeprecationWarning,
    stacklevel=2,
)
