# This file is dual licensed under the terms of the Apache License, Version
# 2.0, and the BSD License. See the LICENSE file in the root of this repository
# for complete details.

from __future__ import annotations

import typing

from cryptography.hazmat.primitives.asymmetric import (
    dh,
    dsa,
    ec,
    ed448,
    ed25519,
    mldsa,
    mlkem,
    rsa,
    x448,
    x25519,
)

# Every asymmetric key type
PublicKeyTypes = typing.Union[
    dh.DHPublicKey,
    dsa.DSAPublicKey,
    rsa.RSAPublicKey,
    ec.EllipticCurvePublicKey,
    ed25519.Ed25519PublicKey,
    ed448.Ed448PublicKey,
    mldsa.MLDSA44PublicKey,
    mldsa.MLDSA65PublicKey,
    mldsa.MLDSA87PublicKey,
    mlkem.MLKEM768PublicKey,
    mlkem.MLKEM1024PublicKey,
    x25519.X25519PublicKey,
    x448.X448PublicKey,
]
# Every asymmetric key type
PrivateKeyTypes = typing.Union[
    dh.DHPrivateKey,
    ed25519.Ed25519PrivateKey,
    ed448.Ed448PrivateKey,
    mldsa.MLDSA44PrivateKey,
    mldsa.MLDSA65PrivateKey,
    mldsa.MLDSA87PrivateKey,
    mlkem.MLKEM768PrivateKey,
    mlkem.MLKEM1024PrivateKey,
    rsa.RSAPrivateKey,
    dsa.DSAPrivateKey,
    ec.EllipticCurvePrivateKey,
    x25519.X25519PrivateKey,
    x448.X448PrivateKey,
]
# Just the key types we allow to be used for x509 signing. This mirrors
# the certificate public key types
CertificateIssuerPrivateKeyTypes = typing.Union[
    ed25519.Ed25519PrivateKey,
    ed448.Ed448PrivateKey,
    rsa.RSAPrivateKey,
    dsa.DSAPrivateKey,
    ec.EllipticCurvePrivateKey,
    mldsa.MLDSA44PrivateKey,
    mldsa.MLDSA65PrivateKey,
    mldsa.MLDSA87PrivateKey,
]
# Just the key types we allow to be used for x509 signing. This mirrors
# the certificate private key types
CertificateIssuerPublicKeyTypes = typing.Union[
    dsa.DSAPublicKey,
    rsa.RSAPublicKey,
    ec.EllipticCurvePublicKey,
    ed25519.Ed25519PublicKey,
    ed448.Ed448PublicKey,
    mldsa.MLDSA44PublicKey,
    mldsa.MLDSA65PublicKey,
    mldsa.MLDSA87PublicKey,
]
# This type removes DHPublicKey. x448/x25519 can be a public key
# but cannot be used in signing so they are allowed here.
CertificatePublicKeyTypes = typing.Union[
    dsa.DSAPublicKey,
    rsa.RSAPublicKey,
    ec.EllipticCurvePublicKey,
    ed25519.Ed25519PublicKey,
    ed448.Ed448PublicKey,
    mldsa.MLDSA44PublicKey,
    mldsa.MLDSA65PublicKey,
    mldsa.MLDSA87PublicKey,
    x25519.X25519PublicKey,
    x448.X448PublicKey,
]
