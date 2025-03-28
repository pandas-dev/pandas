from typing import Any, Dict, List, Optional

from moto.core.utils import iso_8601_datetime_with_milliseconds, utcnow

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
        fingerprint: str,
        material: Optional[str],
        material_public: str,
        tags: Dict[str, str],
        ec2_backend: Any,
    ):
        self.id = random_key_pair_id()
        self.name = name
        self.fingerprint = fingerprint  # public key fingerprint
        self.material = material  # PEM encoded private key
        self.material_public = material_public  # public key in OpenSSH format
        self.create_time = utcnow()
        self.ec2_backend = ec2_backend
        self.add_tags(tags or {})

    @property
    def created_iso_8601(self) -> str:
        return iso_8601_datetime_with_milliseconds(self.create_time)

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> str:
        if filter_name == "key-name":
            return self.name
        elif filter_name == "fingerprint":
            return self.fingerprint
        else:
            return super().get_filter_value(filter_name, "DescribeKeyPairs")


class KeyPairBackend:
    def __init__(self) -> None:
        self.keypairs: Dict[str, KeyPair] = {}

    def create_key_pair(
        self, name: str, key_type: str, tags: Dict[str, str]
    ) -> KeyPair:
        if name in self.keypairs:
            raise InvalidKeyPairDuplicateError(name)
        if key_type == "ed25519":
            keypair = KeyPair(
                name, **random_ed25519_key_pair(), tags=tags, ec2_backend=self
            )
        else:
            keypair = KeyPair(
                name, **random_rsa_key_pair(), tags=tags, ec2_backend=self
            )

        self.keypairs[name] = keypair
        return keypair

    def delete_key_pair(self, name: str) -> None:
        self.keypairs.pop(name, None)

    def describe_key_pairs(
        self, key_names: List[str], filters: Any = None
    ) -> List[KeyPair]:
        if any(key_names):
            results = [
                keypair
                for keypair in self.keypairs.values()
                if keypair.name in key_names
            ]
            if len(key_names) > len(results):
                unknown_keys = set(key_names) - set(results)  # type: ignore
                raise InvalidKeyPairNameError(unknown_keys)
        else:
            results = list(self.keypairs.values())

        if filters:
            return generic_filter(filters, results)
        else:
            return results

    def import_key_pair(
        self, key_name: str, public_key_material: str, tags: Dict[str, str]
    ) -> KeyPair:
        if key_name in self.keypairs:
            raise InvalidKeyPairDuplicateError(key_name)

        try:
            public_key = public_key_parse(public_key_material)
        except ValueError:
            raise InvalidKeyPairFormatError()

        hash_constructor = select_hash_algorithm(public_key)
        fingerprint = public_key_fingerprint(public_key, hash_constructor)
        keypair = KeyPair(
            key_name,
            material_public=public_key_material,
            material=None,
            fingerprint=fingerprint,
            tags=tags,
            ec2_backend=self,
        )
        self.keypairs[key_name] = keypair
        return keypair
