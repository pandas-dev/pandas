from typing import Any, Optional

from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PublicKey

from moto.core.utils import utcnow

from ..exceptions import (
    InvalidKeyPairDuplicateError,
    InvalidKeyPairFormatError,
    InvalidKeyPairNameError,
)
from ..utils import (
    generic_filter,
    public_key_fingerprint,
    public_key_parse,
    random_ed25519_key_pair,
    random_key_pair_id,
    random_rsa_key_pair,
    select_hash_algorithm,
)
from .core import TaggedEC2Resource


class KeyPair(TaggedEC2Resource):
    def __init__(
        self,
        name: str,
        key_type: str,
        fingerprint: str,
        material: bytes,
        material_public: bytes,
        tags: dict[str, str],
        ec2_backend: Any,
    ):
        self.id = random_key_pair_id()
        self.key_name = name
        self.key_type = key_type
        self.key_fingerprint = fingerprint  # public key fingerprint
        self.material = material  # PEM encoded private key
        self.material_public = material_public  # public key in OpenSSH format
        self.create_time = utcnow()
        self.ec2_backend = ec2_backend
        self.add_tags(tags or {})

    @property
    def key_material(self) -> str:
        return self.material.decode("utf-8")

    @property
    def public_key(self) -> str:
        return self.material_public.decode("utf-8")

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> str:
        if filter_name == "key-name":
            return self.key_name
        elif filter_name == "key-pair-id":
            return self.id
        elif filter_name == "fingerprint":
            return self.key_fingerprint
        else:
            return super().get_filter_value(filter_name, "DescribeKeyPairs")


class KeyPairBackend:
    def __init__(self) -> None:
        self.key_pairs: dict[str, KeyPair] = {}

    def create_key_pair(
        self, name: str, key_type: str, tags: dict[str, str]
    ) -> KeyPair:
        if name in self.key_pairs:
            raise InvalidKeyPairDuplicateError(name)
        if key_type == "ed25519":
            key_details = random_ed25519_key_pair()
        else:
            key_details = random_rsa_key_pair()
        keypair = KeyPair(name, key_type, **key_details, tags=tags, ec2_backend=self)
        self.key_pairs[name] = keypair
        return keypair

    def delete_key_pair(self, name: str) -> Optional[KeyPair]:
        return self.key_pairs.pop(name, None)

    def describe_key_pairs(
        self, key_names: list[str], filters: Any = None
    ) -> list[KeyPair]:
        if any(key_names):
            results = [
                keypair
                for keypair in self.key_pairs.values()
                if keypair.key_name in key_names
            ]
            if len(key_names) > len(results):
                unknown_keys = set(key_names) - set(results)  # type: ignore
                raise InvalidKeyPairNameError(unknown_keys)
        else:
            results = list(self.key_pairs.values())

        if filters:
            return generic_filter(filters, results)
        else:
            return results

    def import_key_pair(
        self, key_name: str, public_key_material: bytes, tags: dict[str, str]
    ) -> KeyPair:
        if key_name in self.key_pairs:
            raise InvalidKeyPairDuplicateError(key_name)

        try:
            public_key = public_key_parse(public_key_material)
        except ValueError:
            raise InvalidKeyPairFormatError()

        if isinstance(public_key, Ed25519PublicKey):
            key_type = "ed25519"
        else:
            key_type = "rsa"

        hash_constructor = select_hash_algorithm(public_key)
        fingerprint = public_key_fingerprint(public_key, hash_constructor)
        keypair = KeyPair(
            key_name,
            key_type,
            material_public=public_key_material,
            material=b"",
            fingerprint=fingerprint,
            tags=tags,
            ec2_backend=self,
        )
        self.key_pairs[key_name] = keypair
        return keypair
