# Copyright 2017 Google LLC
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

"""
RSA cryptography signer and verifier.

This file provides a shared wrapper, that defers to _python_rsa or _cryptography_rsa
for implmentations using different third party libraries
"""

from cryptography.hazmat.primitives.asymmetric.rsa import RSAPrivateKey
from cryptography.hazmat.primitives.asymmetric.rsa import RSAPublicKey

from google.auth import _helpers
from google.auth.crypt import _cryptography_rsa
from google.auth.crypt import base

RSA_KEY_MODULE_PREFIX = "rsa.key"


class RSAVerifier(base.Verifier):
    """Verifies RSA cryptographic signatures using public keys.

    Args:
        public_key (Union["rsa.key.PublicKey", cryptography.hazmat.primitives.asymmetric.rsa.RSAPublicKey]):
            The public key used to verify signatures.
    Raises:
        ImportError: if called with an rsa.key.PublicKey, when the rsa library is not installed
        ValueError: if an unrecognized public key is provided
    """

    def __init__(self, public_key):
        module_str = public_key.__class__.__module__
        if isinstance(public_key, RSAPublicKey):
            impl_lib = _cryptography_rsa
        elif module_str.startswith(RSA_KEY_MODULE_PREFIX):
            from google.auth.crypt import _python_rsa

            impl_lib = _python_rsa
        else:
            raise ValueError(f"unrecognized public key type: {type(public_key)}")
        self._impl = impl_lib.RSAVerifier(public_key)

    @_helpers.copy_docstring(base.Verifier)
    def verify(self, message, signature):
        return self._impl.verify(message, signature)

    @classmethod
    def from_string(cls, public_key):
        """Construct a Verifier instance from a public key or public
        certificate string.

        Args:
            public_key (Union[str, bytes]): The public key in PEM format or the
                x509 public key certificate.

        Returns:
            google.auth.crypt.Verifier: The constructed verifier.

        Raises:
            ValueError: If the public_key can't be parsed.
        """
        instance = cls.__new__(cls)
        instance._impl = _cryptography_rsa.RSAVerifier.from_string(public_key)
        return instance


class RSASigner(base.Signer, base.FromServiceAccountMixin):
    """Signs messages with an RSA private key.

    Args:
        private_key (Union["rsa.key.PrivateKey", cryptography.hazmat.primitives.asymmetric.rsa.RSAPrivateKey]):
            The private key to sign with.
        key_id (str): Optional key ID used to identify this private key. This
            can be useful to associate the private key with its associated
            public key or certificate.

    Raises:
        ImportError: if called with an rsa.key.PrivateKey, when the rsa library is not installed
        ValueError: if an unrecognized public key is provided
    """

    def __init__(self, private_key, key_id=None):
        module_str = private_key.__class__.__module__
        if isinstance(private_key, RSAPrivateKey):
            impl_lib = _cryptography_rsa
        elif module_str.startswith(RSA_KEY_MODULE_PREFIX):
            from google.auth.crypt import _python_rsa

            impl_lib = _python_rsa
        else:
            raise ValueError(f"unrecognized private key type: {type(private_key)}")
        self._impl = impl_lib.RSASigner(private_key, key_id=key_id)

    @property  # type: ignore
    @_helpers.copy_docstring(base.Signer)
    def key_id(self):
        return self._impl.key_id

    @_helpers.copy_docstring(base.Signer)
    def sign(self, message):
        return self._impl.sign(message)

    @classmethod
    def from_string(cls, key, key_id=None):
        """Construct a Signer instance from a private key in PEM format.

        Args:
            key (str): Private key in PEM format.
            key_id (str): An optional key id used to identify the private key.

        Returns:
            google.auth.crypt.Signer: The constructed signer.

        Raises:
            ValueError: If the key cannot be parsed as PKCS#1 or PKCS#8 in
                PEM format.
        """
        instance = cls.__new__(cls)
        instance._impl = _cryptography_rsa.RSASigner.from_string(key, key_id=key_id)
        return instance
