import warnings

# Import everything from the new stable module
from google.cloud.storage.asyncio.retry.writes_resumption_strategy import *  # noqa

warnings.warn(
    "google.cloud.storage._experimental.asyncio.retry.writes_resumption_strategy has been moved to google.cloud.storage.asyncio.retry.writes_resumption_strategy. "
    "Please update your imports.",
    DeprecationWarning,
    stacklevel=2,
)
