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

"""ECDSA (ES256) verifier and signer that use the ``cryptography`` library.
"""

from google.auth.crypt.es import EsSigner
from google.auth.crypt.es import EsVerifier


class ES256Verifier(EsVerifier):
    """Verifies ECDSA cryptographic signatures using public keys.

    Args:
        public_key (cryptography.hazmat.primitives.asymmetric.ec.ECDSAPublicKey): The public key used to verify
            signatures.
    """

    pass


class ES256Signer(EsSigner):
    """Signs messages with an ECDSA private key.

    Args:
        private_key (
                cryptography.hazmat.primitives.asymmetric.ec.ECDSAPrivateKey):
            The private key to sign with.
        key_id (str): Optional key ID used to identify this private key. This
            can be useful to associate the private key with its associated
            public key or certificate.
    """

    pass
