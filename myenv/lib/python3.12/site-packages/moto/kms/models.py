import json
import os
from collections import defaultdict
from copy import copy
from datetime import datetime, timedelta
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.exceptions import JsonRESTError
from moto.core.utils import unix_time
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import ValidationException
from .utils import (
    RESERVED_ALIASES,
    KeySpec,
    SigningAlgorithm,
    decrypt,
    encrypt,
    generate_key_id,
    generate_master_key,
    generate_private_key,
)


class Grant(BaseModel):
    def __init__(
        self,
        key_id: str,
        name: str,
        grantee_principal: str,
        operations: List[str],
        constraints: Dict[str, Any],
        retiring_principal: str,
    ):
        self.key_id = key_id
        self.name = name
        self.grantee_principal = grantee_principal
        self.retiring_principal = retiring_principal
        self.operations = operations
        self.constraints = constraints
        self.id = mock_random.get_random_hex()
        self.token = mock_random.get_random_hex()

    def to_json(self) -> Dict[str, Any]:
        return {
            "KeyId": self.key_id,
            "GrantId": self.id,
            "Name": self.name,
            "GranteePrincipal": self.grantee_principal,
            "RetiringPrincipal": self.retiring_principal,
            "Operations": self.operations,
            "Constraints": self.constraints,
        }


class Key(CloudFormationModel):
    def __init__(
        self,
        policy: Optional[str],
        key_usage: str,
        key_spec: str,
        description: str,
        account_id: str,
        region: str,
        multi_region: bool = False,
    ):
        self.id = generate_key_id(multi_region)
        self.creation_date = unix_time()
        self.account_id = account_id
        self.region = region
        self.policy = policy or self.generate_default_policy()
        self.key_usage = key_usage
        self.key_state = "Enabled"
        self.description = description or ""
        self.enabled = True
        self.multi_region = multi_region
        self.key_rotation_status = False
        self.deletion_date: Optional[datetime] = None
        self.key_material = generate_master_key()
        self.origin = "AWS_KMS"
        self.key_manager = "CUSTOMER"
        self.key_spec = key_spec or "SYMMETRIC_DEFAULT"
        self.private_key = generate_private_key(self.key_spec)
        self.arn = (
            f"arn:{get_partition(region)}:kms:{region}:{account_id}:key/{self.id}"
        )

        self.grants: Dict[str, Grant] = dict()

    def add_grant(
        self,
        name: str,
        grantee_principal: str,
        operations: List[str],
        constraints: Dict[str, Any],
        retiring_principal: str,
    ) -> Grant:
        grant = Grant(
            self.id,
            name,
            grantee_principal,
            operations,
            constraints=constraints,
            retiring_principal=retiring_principal,
        )
        self.grants[grant.id] = grant
        return grant

    def list_grants(self, grant_id: str) -> List[Grant]:
        grant_ids = [grant_id] if grant_id else self.grants.keys()
        return [grant for _id, grant in self.grants.items() if _id in grant_ids]

    def list_retirable_grants(self, retiring_principal: str) -> List[Grant]:
        return [
            grant
            for grant in self.grants.values()
            if grant.retiring_principal == retiring_principal
        ]

    def revoke_grant(self, grant_id: str) -> None:
        if not self.grants.pop(grant_id, None):
            raise JsonRESTError("NotFoundException", f"Grant ID {grant_id} not found")

    def retire_grant(self, grant_id: str) -> None:
        self.grants.pop(grant_id, None)

    def retire_grant_by_token(self, grant_token: str) -> None:
        self.grants = {
            _id: grant
            for _id, grant in self.grants.items()
            if grant.token != grant_token
        }

    def generate_default_policy(self) -> str:
        return json.dumps(
            {
                "Version": "2012-10-17",
                "Id": "key-default-1",
                "Statement": [
                    {
                        "Sid": "Enable IAM User Permissions",
                        "Effect": "Allow",
                        "Principal": {
                            "AWS": f"arn:{get_partition(self.region)}:iam::{self.account_id}:root"
                        },
                        "Action": "kms:*",
                        "Resource": "*",
                    }
                ],
            }
        )

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @property
    def encryption_algorithms(self) -> Optional[List[str]]:
        if self.key_usage == "SIGN_VERIFY":
            return None
        elif self.key_spec == "SYMMETRIC_DEFAULT":
            return ["SYMMETRIC_DEFAULT"]
        else:
            return ["RSAES_OAEP_SHA_1", "RSAES_OAEP_SHA_256"]

    @property
    def signing_algorithms(self) -> List[str]:
        if self.key_usage == "ENCRYPT_DECRYPT":
            return None  # type: ignore[return-value]
        elif self.key_spec in KeySpec.ecc_key_specs():
            if self.key_spec == KeySpec.ECC_NIST_P384:
                return [SigningAlgorithm.ECDSA_SHA_384.value]
            elif self.key_spec == KeySpec.ECC_NIST_P521:
                return [SigningAlgorithm.ECDSA_SHA_512.value]
            else:
                # key_spec is 'ECC_NIST_P256' or 'ECC_SECG_P256K1'
                return [SigningAlgorithm.ECDSA_SHA_256.value]
        elif self.key_spec in KeySpec.rsa_key_specs():
            return SigningAlgorithm.rsa_signing_algorithms()
        elif self.key_spec == KeySpec.SM2:
            return [SigningAlgorithm.SM2DSA.value]
        else:
            return []

    def to_dict(self) -> Dict[str, Any]:
        key_dict = {
            "KeyMetadata": {
                "AWSAccountId": self.account_id,
                "Arn": self.arn,
                "CreationDate": self.creation_date,
                "CustomerMasterKeySpec": self.key_spec,
                "KeySpec": self.key_spec,
                "Description": self.description,
                "Enabled": self.enabled,
                "EncryptionAlgorithms": self.encryption_algorithms,
                "KeyId": self.id,
                "KeyManager": self.key_manager,
                "KeyUsage": self.key_usage,
                "KeyState": self.key_state,
                "MultiRegion": self.multi_region,
                "Origin": self.origin,
                "SigningAlgorithms": self.signing_algorithms,
            }
        }
        if self.key_state == "PendingDeletion":
            key_dict["KeyMetadata"]["DeletionDate"] = unix_time(self.deletion_date)
        return key_dict

    def delete(self, account_id: str, region_name: str) -> None:
        kms_backends[account_id][region_name].delete_key(self.id)

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-kms-key.html
        return "AWS::KMS::Key"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Key":
        kms_backend = kms_backends[account_id][region_name]
        properties = cloudformation_json["Properties"]

        key = kms_backend.create_key(
            policy=properties["KeyPolicy"],
            key_usage="ENCRYPT_DECRYPT",
            key_spec="SYMMETRIC_DEFAULT",
            description=properties.get("Description"),
            tags=properties.get("Tags", []),
        )
        key.key_rotation_status = properties.get("EnableKeyRotation", False)
        key.enabled = properties.get("Enabled", True)

        return key

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return self.arn
        raise UnformattedGetAttTemplateException()


class KmsBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: Optional[str] = None):
        super().__init__(region_name=region_name, account_id=account_id)  # type: ignore
        self.keys: Dict[str, Key] = {}
        self.key_to_aliases: Dict[str, Set[str]] = defaultdict(set)
        self.tagger = TaggingService(key_name="TagKey", value_name="TagValue")

    def _generate_default_keys(self, alias_name: str) -> Optional[str]:
        """Creates default kms keys"""
        if alias_name in RESERVED_ALIASES:
            key = self.create_key(
                None,
                "ENCRYPT_DECRYPT",
                "SYMMETRIC_DEFAULT",
                "Default key",
                None,
            )
            self.add_alias(key.id, alias_name)
            return key.id
        return None

    def create_key(
        self,
        policy: Optional[str],
        key_usage: str,
        key_spec: str,
        description: str,
        tags: Optional[List[Dict[str, str]]],
        multi_region: bool = False,
    ) -> Key:
        """
        The provided Policy currently does not need to be valid. If it is valid, Moto will perform authorization checks on key-related operations, just like AWS does.

        These authorization checks are quite basic for now. Moto will only throw an AccessDeniedException if the following conditions are met:
         - The principal is set to "*"
         - The resource is set to "*"
         - The Action matches `describe_key`
        """
        if key_spec:
            self.__ensure_valid_key_spec(key_spec)
        key = Key(
            policy,
            key_usage,
            key_spec,
            description,
            self.account_id,
            self.region_name,
            multi_region,
        )
        self.keys[key.id] = key
        if tags is not None and len(tags) > 0:
            self.tag_resource(key.id, tags)
        return key

    # https://docs.aws.amazon.com/kms/latest/developerguide/multi-region-keys-overview.html#mrk-sync-properties
    # In AWS replicas of a key only share some properties with the original key. Some of those properties get updated
    # in all replicas automatically if those properties change in the original key. Also, such properties can not be
    # changed for replicas directly.
    #
    # In our implementation with just create a copy of all the properties once without any protection from change,
    # as the exact implementation is currently infeasible.
    def replicate_key(self, key_id: str, replica_region: str) -> Key:
        # Using copy() instead of deepcopy(), as the latter results in exception:
        #    TypeError: cannot pickle '_cffi_backend.FFI' object
        # Since we only update top level properties, copy() should suffice.
        replica_key = copy(self.keys[key_id])
        replica_key.region = replica_region
        replica_key.arn = replica_key.arn.replace(self.region_name, replica_region)
        to_region_backend = kms_backends[self.account_id][replica_region]
        to_region_backend.keys[replica_key.id] = replica_key
        return replica_key

    def update_key_description(self, key_id: str, description: str) -> None:
        key = self.keys[self.get_key_id(key_id)]
        key.description = description

    def delete_key(self, key_id: str) -> None:
        if key_id in self.keys:
            if key_id in self.key_to_aliases:
                self.key_to_aliases.pop(key_id)
            self.tagger.delete_all_tags_for_resource(key_id)

            self.keys.pop(key_id)

    def describe_key(self, key_id: str) -> Key:
        # allow the different methods (alias, ARN :key/, keyId, ARN alias) to
        # describe key not just KeyId
        key_id = self.get_key_id(key_id)
        if r"alias/" in str(key_id).lower():
            key_id = self.get_key_id_from_alias(key_id)  # type: ignore[assignment]
        return self.keys[self.get_key_id(key_id)]

    def list_keys(self) -> Iterable[Key]:
        return self.keys.values()

    @staticmethod
    def get_key_id(key_id: str) -> str:
        # Allow use of ARN as well as pure KeyId
        if key_id.startswith("arn:") and ":key/" in key_id:
            return key_id.split(":key/")[1]

        return key_id

    @staticmethod
    def get_alias_name(alias_name: str) -> str:
        # Allow use of ARN as well as alias name
        if alias_name.startswith("arn:") and ":alias/" in alias_name:
            return "alias/" + alias_name.split(":alias/")[1]

        return alias_name

    def any_id_to_key_id(self, key_id: str) -> str:
        """Go from any valid key ID to the raw key ID.

        Acceptable inputs:
        - raw key ID
        - key ARN
        - alias name
        - alias ARN
        """
        key_id = self.get_alias_name(key_id)
        key_id = self.get_key_id(key_id)
        if key_id.startswith("alias/"):
            key_id = self.get_key_id(self.get_key_id_from_alias(key_id))  # type: ignore[arg-type]
        return key_id

    def alias_exists(self, alias_name: str) -> bool:
        for aliases in self.key_to_aliases.values():
            if alias_name in aliases:
                return True

        return False

    def add_alias(self, target_key_id: str, alias_name: str) -> None:
        raw_key_id = self.get_key_id(target_key_id)
        self.key_to_aliases[raw_key_id].add(alias_name)

    def delete_alias(self, alias_name: str) -> None:
        """Delete the alias."""
        for aliases in self.key_to_aliases.values():
            if alias_name in aliases:
                aliases.remove(alias_name)

    def get_all_aliases(self) -> Dict[str, Set[str]]:
        return self.key_to_aliases

    def get_key_id_from_alias(self, alias_name: str) -> Optional[str]:
        for key_id, aliases in dict(self.key_to_aliases).items():
            if alias_name in ",".join(aliases):
                return key_id
        if alias_name in RESERVED_ALIASES:
            return self._generate_default_keys(alias_name)
        return None

    def enable_key_rotation(self, key_id: str) -> None:
        self.keys[self.get_key_id(key_id)].key_rotation_status = True

    def disable_key_rotation(self, key_id: str) -> None:
        self.keys[self.get_key_id(key_id)].key_rotation_status = False

    def get_key_rotation_status(self, key_id: str) -> bool:
        return self.keys[self.get_key_id(key_id)].key_rotation_status

    def put_key_policy(self, key_id: str, policy: str) -> None:
        self.keys[self.get_key_id(key_id)].policy = policy

    def get_key_policy(self, key_id: str) -> str:
        return self.keys[self.get_key_id(key_id)].policy

    def disable_key(self, key_id: str) -> None:
        self.keys[key_id].enabled = False
        self.keys[key_id].key_state = "Disabled"

    def enable_key(self, key_id: str) -> None:
        self.keys[key_id].enabled = True
        self.keys[key_id].key_state = "Enabled"

    def cancel_key_deletion(self, key_id: str) -> None:
        self.keys[key_id].key_state = "Disabled"
        self.keys[key_id].deletion_date = None

    def schedule_key_deletion(self, key_id: str, pending_window_in_days: int) -> float:  # type: ignore[return]
        if 7 <= pending_window_in_days <= 30:
            self.keys[key_id].enabled = False
            self.keys[key_id].key_state = "PendingDeletion"
            self.keys[key_id].deletion_date = datetime.now() + timedelta(
                days=pending_window_in_days
            )
            return unix_time(self.keys[key_id].deletion_date)

    def encrypt(
        self, key_id: str, plaintext: bytes, encryption_context: Dict[str, str]
    ) -> Tuple[bytes, str]:
        key_id = self.any_id_to_key_id(key_id)

        ciphertext_blob = encrypt(
            master_keys=self.keys,
            key_id=key_id,
            plaintext=plaintext,
            encryption_context=encryption_context,
        )
        arn = self.keys[key_id].arn
        return ciphertext_blob, arn

    def decrypt(
        self, ciphertext_blob: bytes, encryption_context: Dict[str, str]
    ) -> Tuple[bytes, str]:
        plaintext, key_id = decrypt(
            master_keys=self.keys,
            ciphertext_blob=ciphertext_blob,
            encryption_context=encryption_context,
        )
        arn = self.keys[key_id].arn
        return plaintext, arn

    def re_encrypt(
        self,
        ciphertext_blob: bytes,
        source_encryption_context: Dict[str, str],
        destination_key_id: str,
        destination_encryption_context: Dict[str, str],
    ) -> Tuple[bytes, str, str]:
        destination_key_id = self.any_id_to_key_id(destination_key_id)

        plaintext, decrypting_arn = self.decrypt(
            ciphertext_blob=ciphertext_blob,
            encryption_context=source_encryption_context,
        )
        new_ciphertext_blob, encrypting_arn = self.encrypt(
            key_id=destination_key_id,
            plaintext=plaintext,
            encryption_context=destination_encryption_context,
        )
        return new_ciphertext_blob, decrypting_arn, encrypting_arn

    def generate_data_key(
        self,
        key_id: str,
        encryption_context: Dict[str, str],
        number_of_bytes: int,
        key_spec: str,
    ) -> Tuple[bytes, bytes, str]:
        key_id = self.any_id_to_key_id(key_id)

        if key_spec:
            # Note: Actual validation of key_spec is done in kms.responses
            if key_spec == "AES_128":
                plaintext_len = 16
            else:
                plaintext_len = 32
        else:
            plaintext_len = number_of_bytes

        plaintext = os.urandom(plaintext_len)

        ciphertext_blob, arn = self.encrypt(
            key_id=key_id, plaintext=plaintext, encryption_context=encryption_context
        )

        return plaintext, ciphertext_blob, arn

    def list_resource_tags(self, key_id_or_arn: str) -> Dict[str, List[Dict[str, str]]]:
        key_id = self.get_key_id(key_id_or_arn)
        if key_id in self.keys:
            return self.tagger.list_tags_for_resource(key_id)
        raise JsonRESTError(
            "NotFoundException",
            "The request was rejected because the specified entity or resource could not be found.",
        )

    def tag_resource(self, key_id_or_arn: str, tags: List[Dict[str, str]]) -> None:
        key_id = self.get_key_id(key_id_or_arn)
        if key_id in self.keys:
            self.tagger.tag_resource(key_id, tags)
            return
        raise JsonRESTError(
            "NotFoundException",
            "The request was rejected because the specified entity or resource could not be found.",
        )

    def untag_resource(self, key_id_or_arn: str, tag_names: List[str]) -> None:
        key_id = self.get_key_id(key_id_or_arn)
        if key_id in self.keys:
            self.tagger.untag_resource_using_names(key_id, tag_names)
            return
        raise JsonRESTError(
            "NotFoundException",
            "The request was rejected because the specified entity or resource could not be found.",
        )

    def create_grant(
        self,
        key_id: str,
        grantee_principal: str,
        operations: List[str],
        name: str,
        constraints: Dict[str, Any],
        retiring_principal: str,
    ) -> Tuple[str, str]:
        key = self.describe_key(key_id)
        grant = key.add_grant(
            name,
            grantee_principal,
            operations,
            constraints=constraints,
            retiring_principal=retiring_principal,
        )
        return grant.id, grant.token

    def list_grants(self, key_id: str, grant_id: str) -> List[Grant]:
        key = self.describe_key(key_id)
        return key.list_grants(grant_id)

    def list_retirable_grants(self, retiring_principal: str) -> List[Grant]:
        grants = []
        for key in self.keys.values():
            grants.extend(key.list_retirable_grants(retiring_principal))
        return grants

    def revoke_grant(self, key_id: str, grant_id: str) -> None:
        key = self.describe_key(key_id)
        key.revoke_grant(grant_id)

    def retire_grant(self, key_id: str, grant_id: str, grant_token: str) -> None:
        if grant_token:
            for key in self.keys.values():
                key.retire_grant_by_token(grant_token)
        else:
            key = self.describe_key(key_id)
            key.retire_grant(grant_id)

    def __ensure_valid_sign_and_verify_key(self, key: Key) -> None:
        if key.key_usage != "SIGN_VERIFY":
            raise ValidationException(
                (
                    "1 validation error detected: Value '{key_id}' at 'KeyId' failed "
                    "to satisfy constraint: Member must point to a key with usage: 'SIGN_VERIFY'"
                ).format(key_id=key.id)
            )

    def __ensure_valid_signing_algorithm(
        self, key: Key, signing_algorithm: str
    ) -> None:
        if signing_algorithm not in key.signing_algorithms:
            raise ValidationException(
                (
                    "1 validation error detected: Value '{signing_algorithm}' at 'SigningAlgorithm' failed "
                    "to satisfy constraint: Member must satisfy enum value set: "
                    "{valid_sign_algorithms}"
                ).format(
                    signing_algorithm=signing_algorithm,
                    valid_sign_algorithms=key.signing_algorithms,
                )
            )

    def __ensure_valid_key_spec(self, key_spec: str) -> None:
        if key_spec not in KeySpec.key_specs():
            raise ValidationException(
                (
                    "1 validation error detected: Value '{key_spec}' at 'KeySpec' failed "
                    "to satisfy constraint: Member must satisfy enum value set: "
                    "{valid_key_specs}"
                ).format(key_spec=key_spec, valid_key_specs=KeySpec.key_specs())
            )

    def sign(
        self, key_id: str, message: bytes, signing_algorithm: str
    ) -> Tuple[str, bytes, str]:
        """
        Sign message using generated private key.

        - grant_tokens are not implemented
        """
        key = self.describe_key(key_id)

        self.__ensure_valid_sign_and_verify_key(key)
        self.__ensure_valid_signing_algorithm(key, signing_algorithm)

        signature = key.private_key.sign(message, signing_algorithm)

        return key.arn, signature, signing_algorithm

    def verify(
        self, key_id: str, message: bytes, signature: bytes, signing_algorithm: str
    ) -> Tuple[str, bool, str]:
        """
        Verify message using public key from generated private key.

        - grant_tokens are not implemented
        - The MessageType-parameter DIGEST is not yet implemented
        """
        key = self.describe_key(key_id)

        self.__ensure_valid_sign_and_verify_key(key)
        self.__ensure_valid_signing_algorithm(key, signing_algorithm)

        if signing_algorithm not in key.signing_algorithms:
            raise ValidationException(
                (
                    "1 validation error detected: Value '{signing_algorithm}' at 'SigningAlgorithm' failed "
                    "to satisfy constraint: Member must satisfy enum value set: "
                    "{valid_sign_algorithms}"
                ).format(
                    signing_algorithm=signing_algorithm,
                    valid_sign_algorithms=key.signing_algorithms,
                )
            )

        return (
            key.arn,
            key.private_key.verify(message, signature, signing_algorithm),
            signing_algorithm,
        )

    def get_public_key(self, key_id: str) -> Tuple[Key, bytes]:
        key = self.describe_key(key_id)
        return key, key.private_key.public_key()


kms_backends = BackendDict(KmsBackend, "kms")
