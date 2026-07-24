from .base import DIRKey as DIRKey
from .cryptography_backend import (
    CryptographyAESKey as CryptographyAESKey,
    CryptographyECKey as CryptographyECKey,
    CryptographyHMACKey as CryptographyHMACKey,
    CryptographyRSAKey as CryptographyRSAKey,
)
from .ecdsa_backend import ECDSAECKey as ECDSAECKey
from .native import HMACKey as NativeHMACKey, get_random_bytes as get_random_bytes
from .rsa_backend import RSAKey as BackendRSAKey

# python-jose relies on importing from cryptography_backend
# then falling back on other imports
# these are all the potential options
AESKey: type[CryptographyAESKey] | None
HMACKey: type[CryptographyHMACKey | NativeHMACKey]
RSAKey: type[CryptographyRSAKey | BackendRSAKey] | None
ECKey: type[CryptographyECKey | ECDSAECKey]
