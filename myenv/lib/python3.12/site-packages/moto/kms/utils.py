import io
import os
import struct
from abc import ABCMeta, abstractmethod
from collections import namedtuple
from enum import Enum
from typing import Any, Dict, List, Tuple

from cryptography.exceptions import InvalidSignature
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives._asymmetric import AsymmetricPadding
from cryptography.hazmat.primitives.asymmetric import ec, padding, rsa
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes

from moto.moto_api._internal import mock_random

from .exceptions import (
    AccessDeniedException,
    InvalidCiphertextException,
    NotFoundException,
    ValidationException,
)

MASTER_KEY_LEN = 32
KEY_ID_LEN = 36
IV_LEN = 12
TAG_LEN = 16
HEADER_LEN = KEY_ID_LEN + IV_LEN + TAG_LEN
# NOTE: This is just a simple binary format. It is not what KMS actually does.
CIPHERTEXT_HEADER_FORMAT = f">{KEY_ID_LEN}s{IV_LEN}s{TAG_LEN}s"
Ciphertext = namedtuple("Ciphertext", ("key_id", "iv", "ciphertext", "tag"))

RESERVED_ALIASE_TARGET_KEY_IDS = {
    # NOTE: These would technically differ across account, but in that they are
    # out of customer control, testing that they are different would be redundant.
    "alias/aws/acm": "4f58743d-e279-4214-9270-8cc28277958d",
    "alias/aws/dynamodb": "7e6aa0ea-15a4-4e72-8b32-58e46e776888",
    "alias/aws/ebs": "7adeb491-68c9-4a5b-86ec-a86ce5364094",
    "alias/aws/elasticfilesystem": "0ef0f111-cdc8-4dda-b0bc-bf625bd5f154",
    "alias/aws/es": "3c7c1880-c353-4cea-9866-d8bc12f05573",
    "alias/aws/glue": "90fd783f-e582-4cc2-a207-672ee67f8d58",
    "alias/aws/kinesisvideo": "7fd4bff3-6eb7-4283-8f11-a7e0a793a181",
    "alias/aws/lambda": "ff9c4f27-2f29-4d9b-bf38-02f88b52a70c",
    "alias/aws/rds": "f5f30938-abed-41a2-a0f6-5482d02a2489",
    "alias/aws/redshift": "dcdae9aa-593a-4e0b-9153-37325591901f",
    "alias/aws/s3": "8c3faf07-f43c-4d11-abdb-9183079214c7",
    "alias/aws/secretsmanager": "fee5173a-3972-428e-ae73-cd4c2a371222",
    "alias/aws/ssm": "cb3f6250-5078-48c0-a75f-0290bf47694e",
    "alias/aws/xray": "e9b758eb-6230-4744-93d1-ad3b7d71f2f6",
}

RESERVED_ALIASES = list(RESERVED_ALIASE_TARGET_KEY_IDS.keys())


class KeySpec(str, Enum):
    # Asymmetric key specs
    RSA_2048 = "RSA_2048"
    RSA_3072 = "RSA_3072"
    RSA_4096 = "RSA_4096"
    ECC_NIST_P256 = "ECC_NIST_P256"
    ECC_SECG_P256K1 = "ECC_SECG_P256K1"
    ECC_NIST_P384 = "ECC_NIST_P384"
    ECC_NIST_P521 = "ECC_NIST_P521"
    SM2 = "SM2"  # China Regions only
    # Symmetric key specs
    SYMMETRIC_DEFAULT = "SYMMETRIC_DEFAULT"
    HMAC_224 = "HMAC_224"
    HMAC_256 = "HMAC_256"
    HMAC_284 = "HMAC_384"
    HMAC_512 = "HMAC_512"

    @classmethod
    def key_specs(self) -> List[str]:
        return sorted([item.value for item in KeySpec])

    @classmethod
    def rsa_key_specs(self) -> List[str]:
        return [spec for spec in self.key_specs() if spec.startswith("RSA")]

    @classmethod
    def ecc_key_specs(self) -> List[str]:
        return [spec for spec in self.key_specs() if spec.startswith("ECC")]

    @classmethod
    def hmac_key_specs(self) -> List[str]:
        return [spec for spec in self.key_specs() if spec.startswith("HMAC")]


class SigningAlgorithm(str, Enum):
    # sigingin algorithms for RSA key spec
    RSASSA_PSS_SHA_256 = "RSASSA_PSS_SHA_256"
    RSASSA_PSS_SHA_384 = "RSASSA_PSS_SHA_384"
    RSASSA_PSS_SHA_512 = "RSASSA_PSS_SHA_512"
    RSASSA_PKCS1_V1_5_SHA_256 = "RSASSA_PKCS1_V1_5_SHA_256"
    RSASSA_PKCS1_V1_5_SHA_384 = "RSASSA_PKCS1_V1_5_SHA_384"
    RSASSA_PKCS1_V1_5_SHA_512 = "RSASSA_PKCS1_V1_5_SHA_512"
    # sigining algorithms for ECC_NIST_P256, P256K1 spec
    ECDSA_SHA_256 = "ECDSA_SHA_256"
    # siginging algorithm for ECC_NIST_P384
    ECDSA_SHA_384 = "ECDSA_SHA_384"
    # sigining algorithm for ECC_NIST_P512
    ECDSA_SHA_512 = "ECDSA_SHA_512"
    # sigining algorithm for SM2
    SM2DSA = "SM2DSA"

    @classmethod
    def signing_algorithms(self) -> List[str]:
        return sorted([item.value for item in SigningAlgorithm])

    @classmethod
    def rsa_signing_algorithms(self) -> List[str]:
        return [algo for algo in self.signing_algorithms() if algo.startswith("RSASSA")]

    @classmethod
    def ecc_signing_algorithms(self) -> List[str]:
        return [algo for algo in self.signing_algorithms() if algo.startswith("ECDSA")]


def generate_key_id(multi_region: bool = False) -> str:
    key = str(mock_random.uuid4())
    # https://docs.aws.amazon.com/kms/latest/developerguide/multi-region-keys-overview.html
    # "Notice that multi-Region keys have a distinctive key ID that begins with mrk-. You can use the mrk- prefix to
    # identify MRKs programmatically."
    if multi_region:
        key = "mrk-" + key

    return key


def generate_data_key(number_of_bytes: int) -> bytes:
    """Generate a data key."""
    return os.urandom(number_of_bytes)


def generate_master_key() -> bytes:
    """Generate a master key."""
    return generate_data_key(MASTER_KEY_LEN)


class AbstractPrivateKey(metaclass=ABCMeta):
    @abstractmethod
    def sign(self, message: bytes, signing_algorithm: str) -> bytes:
        raise NotImplementedError

    @abstractmethod
    def verify(self, message: bytes, signature: bytes, signing_algorithm: str) -> bool:
        raise NotImplementedError

    @abstractmethod
    def public_key(self) -> bytes:
        raise NotImplementedError


def validate_signing_algorithm(
    target_algorithm: str, valid_algorithms: List[str]
) -> None:
    if target_algorithm not in valid_algorithms:
        raise ValidationException(
            (
                "1 validation error detected: Value at 'signing_algorithm' failed"
                "to satisfy constraint: Member must satisfy enum value set: {valid_signing_algorithms}"
            ).format(valid_signing_algorithms=valid_algorithms)
        )


def validate_key_spec(target_key_spec: str, valid_key_specs: List[str]) -> None:
    if target_key_spec not in valid_key_specs:
        raise ValidationException(
            (
                "1 validation error detected: Value at 'key_spec' failed "
                "to satisfy constraint: Member must satisfy enum value set: {valid_key_specs}"
            ).format(valid_key_specs=valid_key_specs)
        )


class RSAPrivateKey(AbstractPrivateKey):
    # See https://docs.aws.amazon.com/kms/latest/cryptographic-details/crypto-primitives.html
    __supported_key_sizes = [2048, 3072, 4096]

    def __init__(self, key_size: int):
        if key_size not in self.__supported_key_sizes:
            raise ValidationException(
                (
                    "1 validation error detected: Value at 'key_size' failed "
                    "to satisfy constraint: Member must satisfy enum value set: {supported_key_sizes}"
                ).format(supported_key_sizes=self.__supported_key_sizes)
            )
        self.key_size = key_size
        self.private_key = rsa.generate_private_key(
            public_exponent=65537, key_size=self.key_size
        )

    def __padding_and_hash_algorithm(
        self, signing_algorithm: str
    ) -> Tuple[AsymmetricPadding, hashes.HashAlgorithm]:
        if signing_algorithm == SigningAlgorithm.RSASSA_PSS_SHA_256:
            pad = padding.PSS(
                mgf=padding.MGF1(hashes.SHA256()), salt_length=padding.PSS.MAX_LENGTH
            )  # type: AsymmetricPadding
            algorithm = hashes.SHA256()  # type: Any
        elif signing_algorithm == SigningAlgorithm.RSASSA_PSS_SHA_384:
            pad = padding.PSS(
                mgf=padding.MGF1(hashes.SHA384()), salt_length=padding.PSS.MAX_LENGTH
            )
            algorithm = hashes.SHA384()
        elif signing_algorithm == SigningAlgorithm.RSASSA_PSS_SHA_512:
            pad = padding.PSS(
                mgf=padding.MGF1(hashes.SHA512()), salt_length=padding.PSS.MAX_LENGTH
            )
            algorithm = hashes.SHA512()
        elif signing_algorithm == SigningAlgorithm.RSASSA_PKCS1_V1_5_SHA_256:
            pad = padding.PKCS1v15()
            algorithm = hashes.SHA256()
        elif signing_algorithm == SigningAlgorithm.RSASSA_PKCS1_V1_5_SHA_384:
            pad = padding.PKCS1v15()
            algorithm = hashes.SHA384()
        else:
            pad = padding.PKCS1v15()
            algorithm = hashes.SHA512()
        return pad, algorithm

    def sign(self, message: bytes, signing_algorithm: str) -> bytes:
        validate_signing_algorithm(
            signing_algorithm, SigningAlgorithm.rsa_signing_algorithms()
        )
        pad, hash_algorithm = self.__padding_and_hash_algorithm(signing_algorithm)
        return self.private_key.sign(message, pad, hash_algorithm)

    def verify(self, message: bytes, signature: bytes, signing_algorithm: str) -> bool:
        validate_signing_algorithm(
            signing_algorithm, SigningAlgorithm.rsa_signing_algorithms()
        )
        pad, hash_algorithm = self.__padding_and_hash_algorithm(signing_algorithm)
        public_key = self.private_key.public_key()
        try:
            public_key.verify(signature, message, pad, hash_algorithm)
            return True
        except InvalidSignature:
            return False

    def public_key(self) -> bytes:
        return self.private_key.public_key().public_bytes(
            encoding=serialization.Encoding.DER,
            format=serialization.PublicFormat.SubjectPublicKeyInfo,
        )


class ECDSAPrivateKey(AbstractPrivateKey):
    def __init__(self, key_spec: str):
        validate_key_spec(key_spec, KeySpec.ecc_key_specs())

        if key_spec == KeySpec.ECC_NIST_P256:
            curve = ec.SECP256R1()  # type: ec.EllipticCurve
            valid_signing_algorithms = ["ECDSA_SHA_256"]  # type: List[str]
        elif key_spec == KeySpec.ECC_SECG_P256K1:
            curve = ec.SECP256K1()
            valid_signing_algorithms = ["ECDSA_SHA_256"]
        elif key_spec == KeySpec.ECC_NIST_P384:
            curve = ec.SECP384R1()
            valid_signing_algorithms = ["ECDSA_SHA_384"]
        else:
            curve = ec.SECP521R1()
            valid_signing_algorithms = ["ECDSA_SHA_512"]

        self.private_key = ec.generate_private_key(curve)
        self.valid_signing_algorithms = valid_signing_algorithms

    def __hash_algorithm(self, signing_algorithm: str) -> hashes.HashAlgorithm:
        if signing_algorithm == SigningAlgorithm.ECDSA_SHA_256:
            algorithm = hashes.SHA256()  # type: Any
        elif signing_algorithm == SigningAlgorithm.ECDSA_SHA_384:
            algorithm = hashes.SHA384()
        else:
            algorithm = hashes.SHA512()
        return algorithm

    def sign(self, message: bytes, signing_algorithm: str) -> bytes:
        validate_signing_algorithm(signing_algorithm, self.valid_signing_algorithms)
        hash_algorithm = self.__hash_algorithm(signing_algorithm)
        return self.private_key.sign(message, ec.ECDSA(hash_algorithm))

    def verify(self, message: bytes, signature: bytes, signing_algorithm: str) -> bool:
        validate_signing_algorithm(signing_algorithm, self.valid_signing_algorithms)
        hash_algorithm = self.__hash_algorithm(signing_algorithm)
        public_key = self.private_key.public_key()
        try:
            public_key.verify(signature, message, ec.ECDSA(hash_algorithm))
            return True
        except InvalidSignature:
            return False

    def public_key(self) -> bytes:
        return self.private_key.public_key().public_bytes(
            encoding=serialization.Encoding.DER,
            format=serialization.PublicFormat.SubjectPublicKeyInfo,
        )


def generate_private_key(key_spec: str) -> AbstractPrivateKey:
    """Generate a private key to be used on asymmetric sign/verify."""
    if key_spec == KeySpec.RSA_2048:
        return RSAPrivateKey(key_size=2048)
    elif key_spec == KeySpec.RSA_3072:
        return RSAPrivateKey(key_size=3072)
    elif key_spec == KeySpec.RSA_4096:
        return RSAPrivateKey(key_size=4096)
    elif key_spec in KeySpec.ecc_key_specs():
        return ECDSAPrivateKey(key_spec)
    else:
        return RSAPrivateKey(key_size=2048)


def _serialize_ciphertext_blob(ciphertext: Ciphertext) -> bytes:
    """Serialize Ciphertext object into a ciphertext blob.

    NOTE: This is just a simple binary format. It is not what KMS actually does.
    """
    header = struct.pack(
        CIPHERTEXT_HEADER_FORMAT,
        ciphertext.key_id.encode("utf-8"),
        ciphertext.iv,
        ciphertext.tag,
    )
    return header + ciphertext.ciphertext


def _deserialize_ciphertext_blob(ciphertext_blob: bytes) -> Ciphertext:
    """Deserialize ciphertext blob into a Ciphertext object.

    NOTE: This is just a simple binary format. It is not what KMS actually does.
    """
    header = ciphertext_blob[:HEADER_LEN]
    ciphertext = ciphertext_blob[HEADER_LEN:]
    key_id, iv, tag = struct.unpack(CIPHERTEXT_HEADER_FORMAT, header)
    return Ciphertext(
        key_id=key_id.decode("utf-8"), iv=iv, ciphertext=ciphertext, tag=tag
    )


def _serialize_encryption_context(encryption_context: Dict[str, str]) -> bytes:
    """Serialize encryption context for use a AAD.

    NOTE: This is not necessarily what KMS does, but it retains the same properties.
    """
    aad = io.BytesIO()
    for key, value in sorted(encryption_context.items(), key=lambda x: x[0]):
        aad.write(key.encode("utf-8"))
        aad.write(value.encode("utf-8"))
    return aad.getvalue()


def encrypt(
    master_keys: Dict[str, Any],
    key_id: str,
    plaintext: bytes,
    encryption_context: Dict[str, str],
) -> bytes:
    """Encrypt data using a master key material.

    NOTE: This is not necessarily what KMS does, but it retains the same properties.

    NOTE: This function is NOT compatible with KMS APIs.
    :param dict master_keys: Mapping of a KmsBackend's known master keys
    :param str key_id: Key ID of moto master key
    :param bytes plaintext: Plaintext data to encrypt
    :param dict[str, str] encryption_context: KMS-style encryption context
    :returns: Moto-structured ciphertext blob encrypted under a moto master key in master_keys
    :rtype: bytes
    """
    try:
        key = master_keys[key_id]
    except KeyError:
        is_alias = key_id.startswith("alias/") or ":alias/" in key_id
        id_type = "Alias" if is_alias else "keyId"
        raise NotFoundException(f"{id_type} {key_id} is not found.")

    if plaintext == b"":
        raise ValidationException(
            "1 validation error detected: Value at 'plaintext' failed to satisfy constraint: Member must have length greater than or equal to 1"
        )
    if len(plaintext) > 4096:
        raise ValidationException(
            "1 validation error detected: Value at 'plaintext' failed to satisfy constraint: Member must have length less than or equal to 4096"
        )

    iv = os.urandom(IV_LEN)
    aad = _serialize_encryption_context(encryption_context=encryption_context)

    encryptor = Cipher(
        algorithms.AES(key.key_material), modes.GCM(iv), backend=default_backend()
    ).encryptor()
    encryptor.authenticate_additional_data(aad)
    ciphertext = encryptor.update(plaintext) + encryptor.finalize()
    return _serialize_ciphertext_blob(
        ciphertext=Ciphertext(
            key_id=key_id, iv=iv, ciphertext=ciphertext, tag=encryptor.tag
        )
    )


def decrypt(
    master_keys: Dict[str, Any],
    ciphertext_blob: bytes,
    encryption_context: Dict[str, str],
) -> Tuple[bytes, str]:
    """Decrypt a ciphertext blob using a master key material.

    NOTE: This is not necessarily what KMS does, but it retains the same properties.

    NOTE: This function is NOT compatible with KMS APIs.

    :param dict master_keys: Mapping of a KmsBackend's known master keys
    :param bytes ciphertext_blob: moto-structured ciphertext blob encrypted under a moto master key in master_keys
    :param dict[str, str] encryption_context: KMS-style encryption context
    :returns: plaintext bytes and moto key ID
    :rtype: bytes and str
    """
    try:
        ciphertext = _deserialize_ciphertext_blob(ciphertext_blob=ciphertext_blob)
    except Exception:
        raise InvalidCiphertextException()

    aad = _serialize_encryption_context(encryption_context=encryption_context)

    try:
        key = master_keys[ciphertext.key_id]
    except KeyError:
        raise AccessDeniedException(
            "The ciphertext refers to a customer master key that does not exist, "
            "does not exist in this region, or you are not allowed to access."
        )

    try:
        decryptor = Cipher(
            algorithms.AES(key.key_material),
            modes.GCM(ciphertext.iv, ciphertext.tag),
            backend=default_backend(),
        ).decryptor()
        decryptor.authenticate_additional_data(aad)
        plaintext = decryptor.update(ciphertext.ciphertext) + decryptor.finalize()
    except Exception:
        raise InvalidCiphertextException()

    return plaintext, ciphertext.key_id
