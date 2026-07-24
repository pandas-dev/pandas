from enum import Enum


class ConnectionStateType(str, Enum):
    AVAILABLE = "available"
    DELETED = "deleted"
    DELETING = "deleting"
    DOWN = "down"
    ORDERING = "ordering"
    PENDING = "pending"
    REJECTED = "rejected"
    REQUESTED = "requested"
    UNKNOWN = "unknown"


class LagStateType(str, Enum):
    AVAILABLE = "available"
    DELETED = "deleted"
    DELETING = "deleting"
    DOWN = "down"
    PENDING = "pending"
    REQUESTED = "requested"
    UNKNOWN = "unknown"


class PortEncryptionStatusType(str, Enum):
    UP = "Encryption Up"
    DOWN = "Encryption Down"


class EncryptionModeType(str, Enum):
    NO = "no_encrypt"
    SHOULD = "should_encrypt"
    MUST = "must_encrypt"


class MacSecKeyStateType(str, Enum):
    ASSOCIATING = "associating"
    ASSOCIATED = "associated"
    DISASSOCIATING = "disassociating"
    DISASSOCIATED = "disassociated"
