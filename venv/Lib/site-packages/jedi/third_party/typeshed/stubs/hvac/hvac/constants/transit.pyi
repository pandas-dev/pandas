from collections.abc import Iterable
from re import Pattern

ALLOWED_KEY_TYPES: Iterable[str]
ALLOWED_EXPORT_KEY_TYPES: Iterable[str]
ALLOWED_DATA_KEY_TYPES: Iterable[str]
ALLOWED_DATA_KEY_BITS: Iterable[int]
ALLOWED_HASH_DATA_ALGORITHMS: Iterable[str]
ALLOWED_HASH_DATA_FORMATS: Iterable[str]
ALLOWED_SIGNATURE_ALGORITHMS: Iterable[str]
ALLOWED_MARSHALING_ALGORITHMS: Iterable[str]
ALLOWED_SALT_LENGTHS: Pattern[str]
