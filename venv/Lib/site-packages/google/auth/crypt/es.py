# Copyright 2017 Google Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""ECDSA verifier and signer that use the ``cryptography`` library.
"""

from dataclasses import dataclass
from typing import Any, Dict, Optional, Union

import cryptography.exceptions
from cryptography.hazmat import backends
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import ec
from cryptography.hazmat.primitives.asymmetric import padding
from cryptography.hazmat.primitives.asymmetric.utils import decode_dss_signature
from cryptography.hazmat.primitives.asymmetric.utils import encode_dss_signature
import cryptography.x509

from google.auth import _helpers
from google.auth.crypt import base


_CERTIFICATE_MARKER = b"-----BEGIN CERTIFICATE-----"
_BACKEND = backends.default_backend()
_PADDING = padding.PKCS1v15()


@dataclass
class _ESAttributes:
    """A class that models ECDSA attributes.

    Attributes:
        rs_size (int): Size for ASN.1 r and s size.
        sha_algo (hashes.HashAlgorithm): Hash algorithm.
        algorithm (str): Algorithm name.
    """

    rs_size: int
    sha_algo: hashes.HashAlgorithm
    algorithm: str

    @classmethod
    def from_key(
        cls, key: Union[ec.EllipticCurvePublicKey, ec.EllipticCurvePrivateKey]
    ):
        return cls.from_curve(key.curve)

    @classmethod
    def from_curve(cls, curve: ec.EllipticCurve):
        # ECDSA raw signature has (r||s) format where r,s are two
        # integers of size 32 bytes for P-256 curve and 48 bytes
        # for P-384 curve. For P-256 curve, we use SHA256 hash algo,
        # and for P-384 curve we use SHA384 algo.
        if isinstance(curve, ec.SECP384R1):
            return cls(48, hashes.SHA384(), "ES384")
        else:
            # default to ES256
            return cls(32, hashes.SHA256(), "ES256")


class EsVerifier(base.Verifier):
    """Verifies ECDSA cryptographic signatures using public keys.

    Args:
        public_key (
                cryptography.hazmat.primitives.asymmetric.ec.EllipticCurvePublicKey):
            The public key used to verify signatures.
    """

    def __init__(self, public_key: ec.EllipticCurvePublicKey) -> None:
        self._pubkey = public_key
        self._attributes = _ESAttributes.from_key(public_key)

    @_helpers.copy_docstring(base.Verifier)
    def verify(self, message: bytes, signature: bytes) -> bool:
        # First convert (r||s) raw signature to ASN1 encoded signature.
        sig_bytes = _helpers.to_bytes(signature)
        if len(sig_bytes) != self._attributes.rs_size * 2:
            return False
        r = int.from_bytes(sig_bytes[: self._attributes.rs_size], byteorder="big")
        s = int.from_bytes(sig_bytes[self._attributes.rs_size :], byteorder="big")
        asn1_sig = encode_dss_signature(r, s)

        message = _helpers.to_bytes(message)
        try:
            self._pubkey.verify(asn1_sig, message, ec.ECDSA(self._attributes.sha_algo))
            return True
        except (ValueError, cryptography.exceptions.InvalidSignature):
            return False

    @classmethod
    def from_string(cls, public_key: Union[str, bytes]) -> "EsVerifier":
        """Construct a Verifier instance from a public key or public
        certificate string.

        Args:
            public_key (Union[str, bytes]): The public key in PEM format or the
                x509 public key certificate.

        Returns:
            google.auth.crypt.Verifier: The constructed verifier.

        Raises:
            ValueError: If the public key can't be parsed.
        """
        public_key_data = _helpers.to_bytes(public_key)

        if _CERTIFICATE_MARKER in public_key_data:
            cert = cryptography.x509.load_pem_x509_certificate(
                public_key_data, _BACKEND
            )
            pubkey = cert.public_key()  # type: Any

        else:
            pubkey = serialization.load_pem_public_key(public_key_data, _BACKEND)

        if not isinstance(pubkey, ec.EllipticCurvePublicKey):
            raise TypeError("Expected public key of type EllipticCurvePublicKey")

        return cls(pubkey)


class EsSigner(base.Signer, base.FromServiceAccountMixin):
    """Signs messages with an ECDSA private key.

    Args:
        private_key (
                cryptography.hazmat.primitives.asymmetric.ec.EllipticCurvePrivateKey):
            The private key to sign with.
        key_id (str): Optional key ID used to identify this private key. This
            can be useful to associate the private key with its associated
            public key or certificate.
    """

    def __init__(
        self, private_key: ec.EllipticCurvePrivateKey, key_id: Optional[str] = None
    ) -> None:
        self._key = private_key
        self._key_id = key_id
        self._attributes = _ESAttributes.from_key(private_key)

    @property
    def algorithm(self) -> str:
        """Name of the algorithm used to sign messages.
        Returns:
            str: The algorithm name.
        """
        return self._attributes.algorithm

    @property  # type: ignore
    @_helpers.copy_docstring(base.Signer)
    def key_id(self) -> Optional[str]:
        return self._key_id

    @_helpers.copy_docstring(base.Signer)
    def sign(self, message: bytes) -> bytes:
        message = _helpers.to_bytes(message)
        asn1_signature = self._key.sign(message, ec.ECDSA(self._attributes.sha_algo))

        # Convert ASN1 encoded signature to (r||s) raw signature.
        (r, s) = decode_dss_signature(asn1_signature)
        return r.to_bytes(self._attributes.rs_size, byteorder="big") + s.to_bytes(
            self._attributes.rs_size, byteorder="big"
        )

    @classmethod
    def from_string(
        cls, key: Union[bytes, str], key_id: Optional[str] = None
    ) -> "EsSigner":
        """Construct a RSASigner from a private key in PEM format.

        Args:
            key (Union[bytes, str]): Private key in PEM format.
            key_id (str): An optional key id used to identify the private key.

        Returns:
            google.auth.crypt._cryptography_rsa.RSASigner: The
            constructed signer.

        Raises:
            ValueError: If ``key`` is not ``bytes`` or ``str`` (unicode).
            UnicodeDecodeError: If ``key`` is ``bytes`` but cannot be decoded
                into a UTF-8 ``str``.
            ValueError: If ``cryptography`` "Could not deserialize key data."
        """
        key_bytes = _helpers.to_bytes(key)
        private_key = serialization.load_pem_private_key(
            key_bytes, password=None, backend=_BACKEND
        )

        if not isinstance(private_key, ec.EllipticCurvePrivateKey):
            raise TypeError("Expected private key of type EllipticCurvePrivateKey")

        return cls(private_key, key_id=key_id)

    def __getstate__(self) -> Dict[str, Any]:
        """Pickle helper that serializes the _key attribute."""
        state = self.__dict__.copy()
        state["_key"] = self._key.private_bytes(
            encoding=serialization.Encoding.PEM,
            format=serialization.PrivateFormat.PKCS8,
            encryption_algorithm=serialization.NoEncryption(),
        )
        return state

    def __setstate__(self, state: Dict[str, Any]) -> None:
        """Pickle helper that deserializes the _key attribute."""
        state["_key"] = serialization.load_pem_private_key(state["_key"], None)
        self.__dict__.update(state)
