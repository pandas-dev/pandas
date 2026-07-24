from typing import Literal

from .abstract.attrDef import AttrDef as AttrDef
from .abstract.attribute import (
    Attribute as Attribute,
    OperationalAttribute as OperationalAttribute,
    WritableAttribute as WritableAttribute,
)
from .abstract.cursor import Reader as Reader, Writer as Writer
from .abstract.entry import Entry as Entry, WritableEntry as WritableEntry
from .abstract.objectDef import ObjectDef as ObjectDef
from .core.connection import Connection as Connection
from .core.pooling import ServerPool as ServerPool
from .core.rdns import ReverseDnsSetting as ReverseDnsSetting
from .core.server import Server as Server
from .core.tls import Tls as Tls
from .protocol.rfc4512 import DsaInfo as DsaInfo, SchemaInfo as SchemaInfo
from .utils.config import get_config_parameter as get_config_parameter, set_config_parameter as set_config_parameter
from .version import __description__ as __description__, __status__ as __status__, __url__ as __url__

ANONYMOUS: Literal["ANONYMOUS"]
SIMPLE: Literal["SIMPLE"]
SASL: Literal["SASL"]
NTLM: Literal["NTLM"]

EXTERNAL: Literal["EXTERNAL"]
DIGEST_MD5: Literal["DIGEST-MD5"]
KERBEROS: Literal["GSSAPI"]
GSSAPI: Literal["GSSAPI"]
PLAIN: Literal["PLAIN"]

AUTO_BIND_DEFAULT: Literal["DEFAULT"]
AUTO_BIND_NONE: Literal["NONE"]
AUTO_BIND_NO_TLS: Literal["NO_TLS"]
AUTO_BIND_TLS_BEFORE_BIND: Literal["TLS_BEFORE_BIND"]
AUTO_BIND_TLS_AFTER_BIND: Literal["TLS_AFTER_BIND"]

IP_SYSTEM_DEFAULT: Literal["IP_SYSTEM_DEFAULT"]
IP_V4_ONLY: Literal["IP_V4_ONLY"]
IP_V6_ONLY: Literal["IP_V6_ONLY"]
IP_V4_PREFERRED: Literal["IP_V4_PREFERRED"]
IP_V6_PREFERRED: Literal["IP_V6_PREFERRED"]

BASE: Literal["BASE"]
LEVEL: Literal["LEVEL"]
SUBTREE: Literal["SUBTREE"]

DEREF_NEVER: Literal["NEVER"]
DEREF_SEARCH: Literal["SEARCH"]
DEREF_BASE: Literal["FINDING_BASE"]
DEREF_ALWAYS: Literal["ALWAYS"]

ALL_ATTRIBUTES: Literal["*"]
NO_ATTRIBUTES: Literal["1.1"]
ALL_OPERATIONAL_ATTRIBUTES: Literal["+"]

MODIFY_ADD: Literal["MODIFY_ADD"]
MODIFY_DELETE: Literal["MODIFY_DELETE"]
MODIFY_REPLACE: Literal["MODIFY_REPLACE"]
MODIFY_INCREMENT: Literal["MODIFY_INCREMENT"]

SYNC: Literal["SYNC"]
SAFE_SYNC: Literal["SAFE_SYNC"]
SAFE_RESTARTABLE: Literal["SAFE_RESTARTABLE"]
ASYNC: Literal["ASYNC"]
LDIF: Literal["LDIF"]
RESTARTABLE: Literal["RESTARTABLE"]
REUSABLE: Literal["REUSABLE"]
MOCK_SYNC: Literal["MOCK_SYNC"]
MOCK_ASYNC: Literal["MOCK_ASYNC"]
ASYNC_STREAM: Literal["ASYNC_STREAM"]

NONE: Literal["NO_INFO"]
DSA: Literal["DSA"]
SCHEMA: Literal["SCHEMA"]
ALL: Literal["ALL"]

OFFLINE_EDIR_8_8_8: Literal["EDIR_8_8_8"]
OFFLINE_EDIR_9_1_4: Literal["EDIR_9_1_4"]
OFFLINE_AD_2012_R2: Literal["AD_2012_R2"]
OFFLINE_SLAPD_2_4: Literal["SLAPD_2_4"]
OFFLINE_DS389_1_3_3: Literal["DS389_1_3_3"]

FIRST: Literal["FIRST"]
ROUND_ROBIN: Literal["ROUND_ROBIN"]
RANDOM: Literal["RANDOM"]

HASHED_NONE: Literal["PLAIN"]
HASHED_SHA: Literal["SHA"]
HASHED_SHA256: Literal["SHA256"]
HASHED_SHA384: Literal["SHA384"]
HASHED_SHA512: Literal["SHA512"]
HASHED_MD5: Literal["MD5"]
HASHED_SALTED_SHA: Literal["SALTED_SHA"]
HASHED_SALTED_SHA256: Literal["SALTED_SHA256"]
HASHED_SALTED_SHA384: Literal["SALTED_SHA384"]
HASHED_SALTED_SHA512: Literal["SALTED_SHA512"]
HASHED_SALTED_MD5: Literal["SALTED_MD5"]

NUMERIC_TYPES: tuple[type, ...]
INTEGER_TYPES: tuple[type, ...]
STRING_TYPES: tuple[type, ...]
SEQUENCE_TYPES: tuple[type, ...]
