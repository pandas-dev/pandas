from enum import Enum
from typing import Final

class Domain(Enum):
    GOOGLE_INTERNAL = 1
    PUBLIC = 2

OSS_DOMAIN: Final[Domain]
OSS_MAJOR: Final[int]
OSS_MINOR: Final[int]
OSS_PATCH: Final[int]
OSS_SUFFIX: Final[str]
DOMAIN: Final[Domain]
MAJOR: Final[int]
MINOR: Final[int]
PATCH: Final[int]
SUFFIX: Final[str]

class VersionError(Exception): ...

def ValidateProtobufRuntimeVersion(
    gen_domain: Domain, gen_major: int, gen_minor: int, gen_patch: int, gen_suffix: str, location: str
) -> None: ...
