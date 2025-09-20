import gzip
import hashlib
import json
import pkgutil
from typing import Any, Dict, Iterator, List, MutableMapping, Optional, Tuple, TypeVar
from uuid import UUID

from requests.structures import CaseInsensitiveDict as _CaseInsensitiveDict

DEFAULT_PARTITION = "aws"
REGION_PREFIX_TO_PARTITION = {
    # (region prefix, aws partition)
    "cn-": "aws-cn",
    "us-gov-": "aws-us-gov",
    "us-iso-": "aws-iso",
    "us-isob-": "aws-iso-b",
}
PARTITION_NAMES = list(REGION_PREFIX_TO_PARTITION.values()) + [DEFAULT_PARTITION]
ARN_PARTITION_REGEX = r"^arn:(" + "|".join(sorted(PARTITION_NAMES)) + ")"


def get_partition(region: str) -> str:
    if not region:
        return DEFAULT_PARTITION
    if region in PARTITION_NAMES:
        return region
    for prefix in REGION_PREFIX_TO_PARTITION:
        if region.startswith(prefix):
            return REGION_PREFIX_TO_PARTITION[prefix]
    return DEFAULT_PARTITION


def str2bool(v: Any) -> Optional[bool]:
    if v in ("yes", True, "true", "True", "TRUE", "t", "1"):
        return True
    elif v in ("no", False, "false", "False", "FALSE", "f", "0"):
        return False
    return None


def load_resource(package: str, resource: str) -> Any:
    """
    Open a file, and return the contents as JSON.
    Will try to load a compressed version (resources/file.json.gz) first, if it exists.
    Usage:
    load_resource(__name__, "resources/file.json")
    """
    try:
        compressed_resource = resource
        if not resource.endswith(".gz"):
            compressed_resource = f"{resource}.gz"
        data = gzip.decompress(pkgutil.get_data(package, compressed_resource))  # type: ignore
    except FileNotFoundError:
        data = pkgutil.get_data(package, resource)  # type: ignore
    return json.loads(data)


def load_resource_as_str(package: str, resource: str) -> str:
    return load_resource_as_bytes(package, resource).decode("utf-8")


def load_resource_as_bytes(package: str, resource: str) -> bytes:
    return pkgutil.get_data(package, resource)  # type: ignore


def merge_multiple_dicts(*args: Any) -> Dict[str, Any]:
    result = {}
    for d in args:
        result.update(d)
    return result


def is_valid_uuid(uuid: str, version: int = 4) -> bool:
    try:
        UUID(uuid, version=version)
        return True
    except ValueError:
        return False


RESOURCE_TYPE = TypeVar("RESOURCE_TYPE")


def filter_resources(
    resources: List[RESOURCE_TYPE],
    filters: Any,
    attr_pairs: Tuple[Tuple[str, ...], ...],
) -> List[RESOURCE_TYPE]:
    """
    Used to filter resources. Usually in get and describe apis.
    """
    result = resources.copy()
    for resource in resources:
        for attrs in attr_pairs:
            values = filters.get(attrs[0]) or None
            if values:
                instance = getattr(resource, attrs[1])
                if (len(attrs) <= 2 and instance not in values) or (
                    len(attrs) == 3 and instance.get(attrs[2]) not in values
                ):
                    result.remove(resource)
                    break
            elif attrs[0] in filters:
                # In case the filter exists but the value of the filter is empty, the filter shouldn't match
                result.remove(resource)
                break
    return result


def md5_hash(data: Any = None) -> Any:
    """
    MD5-hashing for non-security usecases.
    Required for Moto to work in FIPS-enabled systems
    """
    args = (data,) if data else ()
    return hashlib.md5(*args, usedforsecurity=False)  # type: ignore


class LowercaseDict(MutableMapping[str, Any]):
    """A dictionary that lowercases all keys"""

    def __init__(self, *args: Any, **kwargs: Any):
        self.store: Dict[str, Any] = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key: str) -> Any:
        return self.store[self._keytransform(key)]

    def __setitem__(self, key: str, value: Any) -> None:
        self.store[self._keytransform(key)] = value

    def __delitem__(self, key: str) -> None:
        del self.store[self._keytransform(key)]

    def __iter__(self) -> Iterator[Any]:
        return iter(self.store)

    def __len__(self) -> int:
        return len(self.store)

    def __repr__(self) -> str:
        return str(self.store)

    def _keytransform(self, key: str) -> str:
        return key.lower()


class CamelToUnderscoresWalker:
    """A class to convert the keys in dict/list hierarchical data structures from CamelCase to snake_case (underscores)"""

    @staticmethod
    def parse(x: Any) -> Any:  # type: ignore[misc]
        if isinstance(x, dict):
            return CamelToUnderscoresWalker.parse_dict(x)
        elif isinstance(x, list):
            return CamelToUnderscoresWalker.parse_list(x)
        else:
            return CamelToUnderscoresWalker.parse_scalar(x)

    @staticmethod
    def parse_dict(x: Dict[str, Any]) -> Dict[str, Any]:  # type: ignore[misc]
        from moto.core.utils import camelcase_to_underscores

        temp = {}
        for key in x.keys():
            temp[camelcase_to_underscores(key)] = CamelToUnderscoresWalker.parse(x[key])
        return temp

    @staticmethod
    def parse_list(x: Any) -> Any:  # type: ignore[misc]
        temp = []
        for i in x:
            temp.append(CamelToUnderscoresWalker.parse(i))
        return temp

    @staticmethod
    def parse_scalar(x: Any) -> Any:  # type: ignore[misc]
        return x


class CaseInsensitiveDict(_CaseInsensitiveDict):  # type: ignore[type-arg]
    """Proxy for requests.structures.CaseInsensitiveDict"""

    pass
