import base64
import binascii
import hashlib
import logging
import re
import sys
from typing import Any, Dict, Iterator, List, Optional, Tuple, Union
from urllib.parse import urlparse

from requests.structures import CaseInsensitiveDict

from moto.settings import S3_IGNORE_SUBDOMAIN_BUCKETNAME

log = logging.getLogger(__name__)


bucket_name_regex = re.compile(r"(.+)\.s3(.*)\.amazonaws.com")
user_settable_fields = {
    "content-md5",
    "content-language",
    "content-type",
    "content-encoding",
    "cache-control",
    "expires",
    "content-disposition",
    "x-robots-tag",
    "x-amz-checksum-algorithm",
    "x-amz-content-sha256",
    "x-amz-content-crc32",
    "x-amz-content-crc32c",
    "x-amz-content-sha1",
    "x-amz-website-redirect-location",
}
ARCHIVE_STORAGE_CLASSES = [
    "GLACIER",
    "DEEP_ARCHIVE",
    "GLACIER_IR",
]
STORAGE_CLASS = [
    "STANDARD",
    "REDUCED_REDUNDANCY",
    "STANDARD_IA",
    "ONEZONE_IA",
    "INTELLIGENT_TIERING",
] + ARCHIVE_STORAGE_CLASSES
LOGGING_SERVICE_PRINCIPAL = "logging.s3.amazonaws.com"


def bucket_name_from_url(url: str) -> Optional[str]:  # type: ignore
    if S3_IGNORE_SUBDOMAIN_BUCKETNAME:
        return None
    domain = urlparse(url).netloc

    if domain.startswith("www."):
        domain = domain[4:]

    if "amazonaws.com" in domain:
        bucket_result = bucket_name_regex.search(domain)
        if bucket_result:
            return bucket_result.groups()[0]
    else:
        if "." in domain:
            return domain.split(".")[0]
        else:
            # No subdomain found.
            return None


# 'owi-common-cf', 'snippets/test.json' = bucket_and_name_from_url('s3://owi-common-cf/snippets/test.json')
def bucket_and_name_from_url(url: str) -> Union[Tuple[str, str], Tuple[None, None]]:
    prefix = "s3://"
    if url.startswith(prefix):
        bucket_name = url[len(prefix) : url.index("/", len(prefix))]
        key = url[url.index("/", len(prefix)) + 1 :]
        return bucket_name, key
    else:
        return None, None


REGION_URL_REGEX = re.compile(
    r"^https?://(s3[-\.](?P<region1>.+)\.amazonaws\.com/(.+)|"
    r"(.+)\.s3[-\.](?P<region2>.+)\.amazonaws\.com)/?"
)


def parse_region_from_url(url: str, use_default_region: bool = True) -> str:
    match = REGION_URL_REGEX.search(url)
    if match:
        region = match.group("region1") or match.group("region2")
    else:
        region = "us-east-1" if use_default_region else None
    return region


def metadata_from_headers(headers: Dict[str, Any]) -> CaseInsensitiveDict:  # type: ignore
    metadata = CaseInsensitiveDict()  # type: ignore
    meta_regex = re.compile(r"^x-amz-meta-([a-zA-Z0-9\-_.]+)$", flags=re.IGNORECASE)
    for header in headers.keys():
        if isinstance(header, str):
            result = meta_regex.match(header)
            meta_key = None
            if result:
                # Check for extra metadata
                meta_key = result.group(0).lower()
            elif header.lower() in user_settable_fields:
                # Check for special metadata that doesn't start with x-amz-meta
                meta_key = header
            if meta_key:
                metadata[meta_key] = (
                    headers[header][0]
                    if isinstance(headers[header], list)
                    else headers[header]
                )
    return metadata


class _VersionedKeyStore(dict):  # type: ignore
    """A simplified/modified version of Django's `MultiValueDict` taken from:
    https://github.com/django/django/blob/70576740b0bb5289873f5a9a9a4e1a26b2c330e5/django/utils/datastructures.py#L282
    """

    def __sgetitem__(self, key: str) -> List[Any]:
        return super().__getitem__(key)

    def pop(self, key: str) -> None:  # type: ignore
        for version in self.getlist(key, []):
            version.dispose()
        super().pop(key)

    def __getitem__(self, key: str) -> Any:
        return self.__sgetitem__(key)[-1]

    def __setitem__(self, key: str, value: Any) -> Any:
        try:
            current = self.__sgetitem__(key)
            current.append(value)
        except (KeyError, IndexError):
            current = [value]

        super().__setitem__(key, current)

    def get(self, key: str, default: Any = None) -> Any:
        try:
            return self[key]
        except (KeyError, IndexError):
            pass
        return default

    def getlist(self, key: str, default: Any = None) -> Any:
        try:
            return self.__sgetitem__(key)
        except (KeyError, IndexError):
            pass
        return default

    def setlist(self, key: Any, list_: Any) -> Any:
        if isinstance(list_, tuple):
            list_ = list(list_)
        elif not isinstance(list_, list):
            list_ = [list_]

        for existing_version in self.getlist(key, []):
            # Dispose of any FakeKeys that we will not keep
            # We should only have FakeKeys here - but we're checking hasattr to be sure
            if existing_version not in list_ and hasattr(existing_version, "dispose"):
                existing_version.dispose()

        super().__setitem__(key, list_)

    def _iteritems(self) -> Iterator[Tuple[str, Any]]:
        for key in self._self_iterable():
            yield key, self[key]

    def _itervalues(self) -> Iterator[Any]:
        for key in self._self_iterable():
            yield self[key]

    def _iterlists(self) -> Iterator[Tuple[str, List[Any]]]:
        for key in self._self_iterable():
            yield key, self.getlist(key)

    def item_size(self) -> int:
        size = 0
        for val in self._self_iterable().values():
            size += sys.getsizeof(val)
        return size

    def _self_iterable(self) -> Dict[str, Any]:
        # to enable concurrency, return a copy, to avoid "dictionary changed size during iteration"
        # TODO: look into replacing with a locking mechanism, potentially
        return dict(self)

    items = iteritems = _iteritems  # type: ignore
    lists = iterlists = _iterlists
    values = itervalues = _itervalues  # type: ignore


def compute_checksum(body: bytes, algorithm: str, encode_base64: bool = True) -> bytes:
    if algorithm == "SHA1":
        hashed_body = _hash(hashlib.sha1, (body,))
    elif algorithm == "CRC32C":
        try:
            import crc32c

            hashed_body = crc32c.crc32c(body).to_bytes(4, "big")
        except:  # noqa: E722 Do not use bare except
            # Optional library Can't be found - just revert to CRC32
            hashed_body = binascii.crc32(body).to_bytes(4, "big")
    elif algorithm == "CRC32":
        hashed_body = binascii.crc32(body).to_bytes(4, "big")
    else:
        hashed_body = _hash(hashlib.sha256, (body,))
    if encode_base64:
        return base64.b64encode(hashed_body)
    else:
        return hashed_body


def _hash(fn: Any, args: Any) -> bytes:
    try:
        return fn(*args, usedforsecurity=False).digest()
    except TypeError:
        # The usedforsecurity-parameter is only available as of Python 3.9
        return fn(*args).digest()


def cors_matches_origin(origin_header: str, allowed_origins: List[str]) -> bool:
    if "*" in allowed_origins:
        return True
    if origin_header in allowed_origins:
        return True
    for allowed in allowed_origins:
        if re.match(allowed.replace(".", "\\.").replace("*", ".*"), origin_header):
            return True
    return False
