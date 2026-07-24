from _typeshed import Unused
from collections.abc import Callable, Sequence
from enum import Enum
from typing import Any, Literal, NamedTuple, TypeVar, overload
from typing_extensions import LiteralString, Self, TypeAlias, deprecated

from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.asymmetric import ec, rsa
from cryptography.hazmat.primitives.asymmetric.ed448 import Ed448PrivateKey as Ed448PrivateKey, Ed448PublicKey as Ed448PublicKey
from cryptography.hazmat.primitives.asymmetric.ed25519 import (
    Ed25519PrivateKey as Ed25519PrivateKey,
    Ed25519PublicKey as Ed25519PublicKey,
)
from cryptography.hazmat.primitives.asymmetric.x448 import X448PrivateKey as X448PrivateKey, X448PublicKey as X448PublicKey
from cryptography.hazmat.primitives.asymmetric.x25519 import (
    X25519PrivateKey as X25519PrivateKey,
    X25519PublicKey as X25519PublicKey,
)
from jwcrypto.common import JWException

_T = TypeVar("_T")

class UnimplementedOKPCurveKey:
    @classmethod
    def generate(cls) -> None: ...
    @classmethod
    def from_public_bytes(cls, *args) -> None: ...
    @classmethod
    def from_private_bytes(cls, *args) -> None: ...

ImplementedOkpCurves: Sequence[str]
priv_bytes: Callable[[bytes], X25519PrivateKey] | None

class _Ed25519_CURVE(NamedTuple):
    pubkey: UnimplementedOKPCurveKey
    privkey: UnimplementedOKPCurveKey

class _Ed448_CURVE(NamedTuple):
    pubkey: UnimplementedOKPCurveKey
    privkey: UnimplementedOKPCurveKey

class _X25519_CURVE(NamedTuple):
    pubkey: UnimplementedOKPCurveKey
    privkey: UnimplementedOKPCurveKey

class _X448_CURVE(NamedTuple):
    pubkey: UnimplementedOKPCurveKey
    privkey: UnimplementedOKPCurveKey

_JWKKeyTypeSupported: TypeAlias = Literal["oct", "RSA", "EC", "OKP"]
JWKTypesRegistry: dict[_JWKKeyTypeSupported, str]

class ParmType(Enum):
    name = "A string with a name"  # pyright: ignore[reportAssignmentType]
    b64 = "Base64url Encoded"
    b64u = "Base64urlUint Encoded"
    unsupported = "Unsupported Parameter"

class JWKParameter(NamedTuple):
    description: str
    public: bool
    required: bool | None
    type: ParmType | None

JWKValuesRegistry: dict[LiteralString, dict[LiteralString, JWKParameter]]
JWKParamsRegistry: dict[LiteralString, JWKParameter]
JWKEllipticCurveRegistry: dict[LiteralString, str]
_JWKUseSupported: TypeAlias = Literal["sig", "enc"]
JWKUseRegistry: dict[_JWKUseSupported, str]
_JWKOperationSupported: TypeAlias = Literal[
    "sign", "verify", "encrypt", "decrypt", "wrapKey", "unwrapKey", "deriveKey", "deriveBits"
]
JWKOperationsRegistry: dict[_JWKOperationSupported, str]
JWKpycaCurveMap: dict[LiteralString, LiteralString]
IANANamedInformationHashAlgorithmRegistry: dict[
    LiteralString,
    hashes.SHA256
    | hashes.SHA384
    | hashes.SHA512
    | hashes.SHA3_224
    | hashes.SHA3_256
    | hashes.SHA3_384
    | hashes.SHA3_512
    | hashes.BLAKE2s
    | hashes.BLAKE2b
    | None,
]

class InvalidJWKType(JWException):
    value: str | None
    def __init__(self, value: str | None = None) -> None: ...

class InvalidJWKUsage(JWException):
    value: str
    use: str
    def __init__(self, use: str, value: str) -> None: ...

class InvalidJWKOperation(JWException):
    op: str
    values: Sequence[str]
    def __init__(self, operation: str, values: Sequence[str]) -> None: ...

class InvalidJWKValue(JWException): ...

class JWK(dict[str, Any]):
    unsafe_skip_rsa_key_validation: bool
    def __init__(self, **kwargs) -> None: ...
    # `kty` and the other keyword arguments are passed as `params` to the called generator
    # function. The possible arguments depend on the value of `kty`.
    # TODO: Add overloads for the individual `kty` values.
    @classmethod
    @overload
    def generate(
        cls,
        *,
        kty: Literal["RSA"],
        public_exponent: int | None = None,
        size: int | None = None,
        kid: str | None = None,
        alg: str | None = None,
        use: _JWKUseSupported | None = None,
        key_ops: list[_JWKOperationSupported] | None = None,
    ) -> Self: ...
    @classmethod
    @overload
    def generate(cls, *, kty: _JWKKeyTypeSupported, **kwargs) -> Self: ...
    def generate_key(self, *, kty: _JWKKeyTypeSupported, **kwargs) -> None: ...
    def import_key(self, **kwargs) -> None: ...
    @classmethod
    def from_json(cls, key) -> Self: ...
    @overload
    def export(self, private_key: bool = True, as_dict: Literal[False] = False) -> str: ...
    @overload
    def export(self, private_key: bool, as_dict: Literal[True]) -> dict[str, Any]: ...
    @overload
    def export(self, *, as_dict: Literal[True]) -> dict[str, Any]: ...
    @overload
    def export_public(self, as_dict: Literal[False] = False) -> str: ...
    @overload
    def export_public(self, as_dict: Literal[True]) -> dict[str, Any]: ...
    @overload
    def export_public(self, as_dict: bool = False) -> str | dict[str, Any]: ...
    @overload
    def export_private(self, as_dict: Literal[False] = False) -> str: ...
    @overload
    def export_private(self, as_dict: Literal[True]) -> dict[str, Any]: ...
    @overload
    def export_private(self, as_dict: bool = False) -> str | dict[str, Any]: ...
    @overload
    def export_symmetric(self, as_dict: Literal[False] = False) -> str: ...
    @overload
    def export_symmetric(self, as_dict: Literal[True]) -> dict[str, Any]: ...
    @overload
    def export_symmetric(self, as_dict: bool = False) -> str | dict[str, Any]: ...
    def public(self) -> Self: ...
    @property
    def has_public(self) -> bool: ...
    @property
    def has_private(self) -> bool: ...
    @property
    def is_symmetric(self) -> bool: ...
    @property
    @deprecated("")
    def key_type(self) -> str | None: ...
    @property
    @deprecated("")
    def key_id(self) -> str | None: ...
    @property
    @deprecated("")
    def key_curve(self) -> str | None: ...
    @deprecated("")
    def get_curve(
        self, arg: str
    ) -> (
        ec.SECP256R1
        | ec.SECP384R1
        | ec.SECP521R1
        | ec.SECP256K1
        | ec.BrainpoolP256R1
        | ec.BrainpoolP384R1
        | ec.BrainpoolP512R1
        | _Ed25519_CURVE
        | _Ed448_CURVE
        | _X25519_CURVE
        | _X448_CURVE
    ): ...
    def get_op_key(
        self, operation: str | None = None, arg: str | None = None
    ) -> str | rsa.RSAPrivateKey | rsa.RSAPublicKey | ec.EllipticCurvePrivateKey | ec.EllipticCurvePublicKey | None: ...
    def import_from_pyca(
        self,
        key: (
            rsa.RSAPrivateKey
            | rsa.RSAPublicKey
            | ec.EllipticCurvePrivateKey
            | ec.EllipticCurvePublicKey
            | Ed25519PrivateKey
            | Ed448PrivateKey
            | X25519PrivateKey
            | Ed25519PublicKey
            | Ed448PublicKey
            | X25519PublicKey
        ),
    ) -> None: ...
    def import_from_pem(self, data: bytes, password: bytes | None = None, kid: str | None = None) -> None: ...
    @overload
    def export_to_pem(self, private_key: Literal[False] = False, password: Unused = False) -> bytes: ...
    @overload
    def export_to_pem(self, private_key: Literal[True], password: bytes | None) -> bytes: ...
    @classmethod
    def from_pyca(
        cls,
        key: (
            rsa.RSAPrivateKey
            | rsa.RSAPublicKey
            | ec.EllipticCurvePrivateKey
            | ec.EllipticCurvePublicKey
            | Ed25519PrivateKey
            | Ed448PrivateKey
            | X25519PrivateKey
            | Ed25519PublicKey
            | Ed448PublicKey
            | X25519PublicKey
        ),
    ) -> Self: ...
    @classmethod
    def from_pem(cls, data: bytes, password: bytes | None = None) -> Self: ...
    def thumbprint(self, hashalg: hashes.HashAlgorithm = ...) -> str: ...
    def thumbprint_uri(self, hname: str = "sha-256") -> str: ...
    @classmethod
    def from_password(cls, password: str) -> Self: ...
    def setdefault(self, key: str, default: _T | None = None) -> _T: ...

class JWKSet(dict[Literal["keys"], set[JWK]]):
    @overload
    def __setitem__(self, key: Literal["keys"], val: JWK) -> None: ...
    @overload
    def __setitem__(self, key: str, val: Any) -> None: ...
    def add(self, elem: JWK) -> None: ...
    @overload
    def export(self, private_keys: bool = True, as_dict: Literal[False] = False) -> str: ...
    @overload
    def export(self, private_keys: bool, as_dict: Literal[True]) -> dict[str, Any]: ...
    @overload
    def export(self, *, as_dict: Literal[True]) -> dict[str, Any]: ...
    def import_keyset(self, keyset: str | bytes) -> None: ...
    @classmethod
    def from_json(cls, keyset: str | bytes) -> Self: ...
    def get_key(self, kid: str) -> JWK | None: ...
    def get_keys(self, kid: str) -> set[JWK]: ...
    def setdefault(self, key: str, default: _T | None = None) -> _T: ...
