import base64
import json
import os
import re
from typing import Any, Dict

from moto.core.responses import BaseResponse
from moto.kms.utils import RESERVED_ALIASE_TARGET_KEY_IDS, RESERVED_ALIASES
from moto.utilities.utils import get_partition

from .exceptions import (
    AlreadyExistsException,
    NotAuthorizedException,
    NotFoundException,
    ValidationException,
)
from .models import KmsBackend, kms_backends
from .policy_validator import validate_policy


class KmsResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="kms")

    def _get_param(self, param_name: str, if_none: Any = None) -> Any:
        params = json.loads(self.body)

        for key in ("Plaintext", "CiphertextBlob", "Message"):
            if key in params:
                params[key] = base64.b64decode(params[key].encode("utf-8"))

        return params.get(param_name, if_none)

    @property
    def kms_backend(self) -> KmsBackend:
        return kms_backends[self.current_account][self.region]

    def _display_arn(self, key_id: str) -> str:
        if key_id.startswith("arn:"):
            return key_id

        if key_id.startswith("alias/"):
            id_type = ""
        else:
            id_type = "key/"

        return f"arn:{get_partition(self.region)}:kms:{self.region}:{self.current_account}:{id_type}{key_id}"

    def _validate_cmk_id(self, key_id: str) -> None:
        """Determine whether a CMK ID exists.

        - raw key ID
        - key ARN
        """
        is_arn = key_id.startswith("arn:") and ":key/" in key_id
        # https://docs.aws.amazon.com/kms/latest/developerguide/multi-region-keys-overview.html
        # "Notice that multi-Region keys have a distinctive key ID that begins with mrk-. You can use the mrk- prefix to
        # identify MRKs programmatically."
        is_raw_key_id = re.match(
            r"^(mrk-)?[A-F0-9]{8}-[A-F0-9]{4}-[A-F0-9]{4}-[A-F0-9]{4}-[A-F0-9]{12}$",
            key_id,
            re.IGNORECASE,
        )

        if not is_arn and not is_raw_key_id:
            raise NotFoundException(f"Invalid keyId {key_id}")

        cmk_id = self.kms_backend.get_key_id(key_id)

        if cmk_id not in self.kms_backend.keys:
            raise NotFoundException(f"Key '{self._display_arn(key_id)}' does not exist")

    def _validate_alias(self, key_id: str) -> None:
        """Determine whether an alias exists.

        - alias name
        - alias ARN
        """
        error = NotFoundException(f"Alias {self._display_arn(key_id)} is not found.")

        is_arn = key_id.startswith("arn:") and ":alias/" in key_id
        is_name = key_id.startswith("alias/")

        if not is_arn and not is_name:
            raise error

        alias_name = self.kms_backend.get_alias_name(key_id)
        cmk_id = self.kms_backend.get_key_id_from_alias(alias_name)
        if cmk_id is None:
            raise error

    def _validate_key_id(self, key_id: str) -> None:
        """Determine whether a key ID exists.

        - raw key ID
        - key ARN
        - alias name
        - alias ARN
        """
        is_alias_arn = key_id.startswith("arn:") and ":alias/" in key_id
        is_alias_name = key_id.startswith("alias/")

        if is_alias_arn or is_alias_name:
            self._validate_alias(key_id)
            return

        self._validate_cmk_id(key_id)

    def _validate_key_policy(self, key_id: str, action: str) -> None:
        """
        Validate whether the specified action is allowed, given the key policy
        """
        key = self.kms_backend.describe_key(self.kms_backend.get_key_id(key_id))
        validate_policy(key, action)

    def create_key(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_CreateKey.html"""
        policy = self._get_param("Policy")
        key_usage = self._get_param("KeyUsage")
        key_spec = self._get_param("KeySpec") or self._get_param(
            "CustomerMasterKeySpec"
        )
        description = self._get_param("Description")
        tags = self._get_param("Tags")
        multi_region = self._get_param("MultiRegion")

        key = self.kms_backend.create_key(
            policy, key_usage, key_spec, description, tags, multi_region
        )
        return json.dumps(key.to_dict())

    def replicate_key(self) -> str:
        key_id = self._get_param("KeyId")
        self._validate_key_id(key_id)
        replica_region = self._get_param("ReplicaRegion")
        replica_key = self.kms_backend.replicate_key(key_id, replica_region)
        return json.dumps(
            {
                "ReplicaKeyMetadata": replica_key.to_dict()["KeyMetadata"],
                "ReplicaPolicy": replica_key.generate_default_policy(),
            }
        )

    def update_key_description(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_UpdateKeyDescription.html"""
        key_id = self._get_param("KeyId")
        description = self._get_param("Description")

        self._validate_cmk_id(key_id)

        self.kms_backend.update_key_description(key_id, description)
        return json.dumps(None)

    def tag_resource(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_TagResource.html"""
        key_id = self._get_param("KeyId")
        tags = self._get_param("Tags")

        self._validate_cmk_id(key_id)

        self.kms_backend.tag_resource(key_id, tags)
        return "{}"

    def untag_resource(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_UntagResource.html"""
        key_id = self._get_param("KeyId")
        tag_names = self._get_param("TagKeys")

        self._validate_cmk_id(key_id)

        self.kms_backend.untag_resource(key_id, tag_names)
        return "{}"

    def list_resource_tags(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_ListResourceTags.html"""
        key_id = self._get_param("KeyId")
        self._validate_cmk_id(key_id)

        tags: Dict[str, Any] = self.kms_backend.list_resource_tags(key_id)
        tags.update({"NextMarker": None, "Truncated": False})
        return json.dumps(tags)

    def describe_key(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_DescribeKey.html"""
        key_id = self._get_param("KeyId")

        self._validate_key_id(key_id)
        self._validate_key_policy(key_id, "kms:DescribeKey")

        key = self.kms_backend.describe_key(self.kms_backend.get_key_id(key_id))

        return json.dumps(key.to_dict())

    def list_keys(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_ListKeys.html"""
        keys = self.kms_backend.list_keys()

        return json.dumps(
            {
                "Keys": [{"KeyArn": key.arn, "KeyId": key.id} for key in keys],
                "NextMarker": None,
                "Truncated": False,
            }
        )

    def create_alias(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_CreateAlias.html"""
        return self._set_alias()

    def update_alias(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_UpdateAlias.html"""
        return self._set_alias(update=True)

    def _set_alias(self, update: bool = False) -> str:
        alias_name = self._get_param("AliasName")
        target_key_id = self._get_param("TargetKeyId")

        if not alias_name.startswith("alias/"):
            raise ValidationException("Invalid identifier")

        if alias_name in RESERVED_ALIASES:
            raise NotAuthorizedException()

        if ":" in alias_name:
            raise ValidationException(
                f"{alias_name} contains invalid characters for an alias"
            )

        if not re.match(r"^[a-zA-Z0-9:/_-]+$", alias_name):
            raise ValidationException(
                f"1 validation error detected: Value '{alias_name}' at 'aliasName' "
                "failed to satisfy constraint: Member must satisfy regular "
                "expression pattern: ^[a-zA-Z0-9:/_-]+$"
            )

        if self.kms_backend.alias_exists(target_key_id):
            raise ValidationException("Aliases must refer to keys. Not aliases")

        if update:
            # delete any existing aliases with that name (should be a no-op if none exist)
            self.kms_backend.delete_alias(alias_name)

        if self.kms_backend.alias_exists(alias_name):
            raise AlreadyExistsException(
                f"An alias with the name arn:aws:kms:{self.region}:{self.current_account}:{alias_name} already exists"
            )

        self._validate_cmk_id(target_key_id)
        self.kms_backend.add_alias(target_key_id, alias_name)

        return json.dumps(None)

    def delete_alias(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_DeleteAlias.html"""
        alias_name = self._get_param("AliasName")

        if not alias_name.startswith("alias/"):
            raise ValidationException("Invalid identifier")

        self._validate_alias(alias_name)

        self.kms_backend.delete_alias(alias_name)

        return json.dumps(None)

    def list_aliases(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_ListAliases.html"""
        region = self.region
        key_id = self._get_param("KeyId")
        if key_id is not None:
            self._validate_key_id(key_id)
            key_id = self.kms_backend.get_key_id(key_id)

        response_aliases = []

        backend_aliases = self.kms_backend.get_all_aliases()
        for target_key_id, aliases in backend_aliases.items():
            for alias_name in aliases:
                # TODO: add creation date and last updated in response_aliases
                response_aliases.append(
                    {
                        "AliasArn": f"arn:{get_partition(region)}:kms:{region}:{self.current_account}:{alias_name}",
                        "AliasName": alias_name,
                        "TargetKeyId": target_key_id,
                    }
                )
        for reserved_alias, target_key_id in RESERVED_ALIASE_TARGET_KEY_IDS.items():
            exsisting = [
                a for a in response_aliases if a["AliasName"] == reserved_alias
            ]
            if not exsisting:
                arn = f"arn:{get_partition(region)}:kms:{region}:{self.current_account}:{reserved_alias}"
                response_aliases.append(
                    {
                        "TargetKeyId": target_key_id,
                        "AliasArn": arn,
                        "AliasName": reserved_alias,
                    }
                )

        if key_id is not None:
            response_aliases = list(
                filter(lambda alias: alias["TargetKeyId"] == key_id, response_aliases)
            )

        return json.dumps({"Truncated": False, "Aliases": response_aliases})

    def create_grant(self) -> str:
        key_id = self._get_param("KeyId")
        grantee_principal = self._get_param("GranteePrincipal")
        retiring_principal = self._get_param("RetiringPrincipal")
        operations = self._get_param("Operations")
        name = self._get_param("Name")
        constraints = self._get_param("Constraints")

        grant_id, grant_token = self.kms_backend.create_grant(
            key_id,
            grantee_principal,
            operations,
            name,
            constraints=constraints,
            retiring_principal=retiring_principal,
        )
        return json.dumps({"GrantId": grant_id, "GrantToken": grant_token})

    def list_grants(self) -> str:
        key_id = self._get_param("KeyId")
        grant_id = self._get_param("GrantId")

        grants = self.kms_backend.list_grants(key_id=key_id, grant_id=grant_id)
        return json.dumps(
            {
                "Grants": [gr.to_json() for gr in grants],
                "GrantCount": len(grants),
                "Truncated": False,
            }
        )

    def list_retirable_grants(self) -> str:
        retiring_principal = self._get_param("RetiringPrincipal")

        grants = self.kms_backend.list_retirable_grants(retiring_principal)
        return json.dumps(
            {
                "Grants": [gr.to_json() for gr in grants],
                "GrantCount": len(grants),
                "Truncated": False,
            }
        )

    def revoke_grant(self) -> str:
        key_id = self._get_param("KeyId")
        grant_id = self._get_param("GrantId")

        self.kms_backend.revoke_grant(key_id, grant_id)
        return "{}"

    def retire_grant(self) -> str:
        key_id = self._get_param("KeyId")
        grant_id = self._get_param("GrantId")
        grant_token = self._get_param("GrantToken")

        self.kms_backend.retire_grant(key_id, grant_id, grant_token)
        return "{}"

    def enable_key_rotation(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_EnableKeyRotation.html"""
        key_id = self._get_param("KeyId")

        self._validate_cmk_id(key_id)

        self.kms_backend.enable_key_rotation(key_id)

        return json.dumps(None)

    def disable_key_rotation(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_EnableKeyRotation.html"""
        key_id = self._get_param("KeyId")

        self._validate_cmk_id(key_id)

        self.kms_backend.disable_key_rotation(key_id)

        return json.dumps(None)

    def get_key_rotation_status(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_GetKeyRotationStatus.html"""
        key_id = self._get_param("KeyId")

        self._validate_cmk_id(key_id)

        rotation_enabled = self.kms_backend.get_key_rotation_status(key_id)

        return json.dumps({"KeyRotationEnabled": rotation_enabled})

    def put_key_policy(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_PutKeyPolicy.html"""
        key_id = self._get_param("KeyId")
        policy_name = self._get_param("PolicyName")
        policy = self._get_param("Policy")
        _assert_default_policy(policy_name)

        self._validate_cmk_id(key_id)

        self.kms_backend.put_key_policy(key_id, policy)

        return json.dumps(None)

    def get_key_policy(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_GetKeyPolicy.html"""
        key_id = self._get_param("KeyId")
        policy_name = self._get_param("PolicyName")
        _assert_default_policy(policy_name)

        self._validate_cmk_id(key_id)

        policy = self.kms_backend.get_key_policy(key_id) or "{}"
        return json.dumps({"Policy": policy})

    def list_key_policies(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_ListKeyPolicies.html"""
        key_id = self._get_param("KeyId")

        self._validate_cmk_id(key_id)

        self.kms_backend.describe_key(key_id)

        return json.dumps({"Truncated": False, "PolicyNames": ["default"]})

    def encrypt(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_Encrypt.html"""
        key_id = self._get_param("KeyId")
        encryption_context = self._get_param("EncryptionContext", {})
        plaintext = self._get_param("Plaintext")

        self._validate_key_id(key_id)

        if isinstance(plaintext, str):
            plaintext = plaintext.encode("utf-8")

        ciphertext_blob, arn = self.kms_backend.encrypt(
            key_id=key_id, plaintext=plaintext, encryption_context=encryption_context
        )
        ciphertext_blob_response = base64.b64encode(ciphertext_blob).decode("utf-8")

        return json.dumps({"CiphertextBlob": ciphertext_blob_response, "KeyId": arn})

    def decrypt(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_Decrypt.html"""
        ciphertext_blob = self._get_param("CiphertextBlob")
        encryption_context = self._get_param("EncryptionContext", {})

        plaintext, arn = self.kms_backend.decrypt(
            ciphertext_blob=ciphertext_blob, encryption_context=encryption_context
        )

        plaintext_response = base64.b64encode(plaintext).decode("utf-8")

        return json.dumps({"Plaintext": plaintext_response, "KeyId": arn})

    def re_encrypt(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_ReEncrypt.html"""
        ciphertext_blob = self._get_param("CiphertextBlob")
        source_encryption_context = self._get_param("SourceEncryptionContext", {})
        destination_key_id = self._get_param("DestinationKeyId")
        destination_encryption_context = self._get_param(
            "DestinationEncryptionContext", {}
        )

        self._validate_key_id(destination_key_id)

        (
            new_ciphertext_blob,
            decrypting_arn,
            encrypting_arn,
        ) = self.kms_backend.re_encrypt(
            ciphertext_blob=ciphertext_blob,
            source_encryption_context=source_encryption_context,
            destination_key_id=destination_key_id,
            destination_encryption_context=destination_encryption_context,
        )

        response_ciphertext_blob = base64.b64encode(new_ciphertext_blob).decode("utf-8")

        return json.dumps(
            {
                "CiphertextBlob": response_ciphertext_blob,
                "KeyId": encrypting_arn,
                "SourceKeyId": decrypting_arn,
            }
        )

    def disable_key(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_DisableKey.html"""
        key_id = self._get_param("KeyId")

        self._validate_cmk_id(key_id)

        self.kms_backend.disable_key(key_id)

        return json.dumps(None)

    def enable_key(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_EnableKey.html"""
        key_id = self._get_param("KeyId")

        self._validate_cmk_id(key_id)

        self.kms_backend.enable_key(key_id)

        return json.dumps(None)

    def cancel_key_deletion(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_CancelKeyDeletion.html"""
        key_id = self._get_param("KeyId")

        self._validate_cmk_id(key_id)

        self.kms_backend.cancel_key_deletion(key_id)

        return json.dumps({"KeyId": key_id})

    def schedule_key_deletion(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_ScheduleKeyDeletion.html"""
        key_id = self._get_param("KeyId")
        if self._get_param("PendingWindowInDays") is None:
            pending_window_in_days = 30
        else:
            pending_window_in_days = self._get_param("PendingWindowInDays")

        self._validate_cmk_id(key_id)

        return json.dumps(
            {
                "KeyId": key_id,
                "DeletionDate": self.kms_backend.schedule_key_deletion(
                    key_id, pending_window_in_days
                ),
            }
        )

    def generate_data_key(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_GenerateDataKey.html"""
        key_id = self._get_param("KeyId")
        encryption_context = self._get_param("EncryptionContext", {})
        number_of_bytes = self._get_param("NumberOfBytes")
        key_spec = self._get_param("KeySpec")

        # Param validation
        self._validate_key_id(key_id)

        if number_of_bytes and (number_of_bytes > 1024 or number_of_bytes < 1):
            raise ValidationException(
                (
                    "1 validation error detected: Value '{number_of_bytes:d}' at 'numberOfBytes' failed "
                    "to satisfy constraint: Member must have value less than or "
                    "equal to 1024"
                ).format(number_of_bytes=number_of_bytes)
            )

        if key_spec and key_spec not in ("AES_256", "AES_128"):
            raise ValidationException(
                (
                    "1 validation error detected: Value '{key_spec}' at 'keySpec' failed "
                    "to satisfy constraint: Member must satisfy enum value set: "
                    "[AES_256, AES_128]"
                ).format(key_spec=key_spec)
            )
        if not key_spec and not number_of_bytes:
            raise ValidationException(
                "Please specify either number of bytes or key spec."
            )

        if key_spec and number_of_bytes:
            raise ValidationException(
                "Please specify either number of bytes or key spec."
            )

        plaintext, ciphertext_blob, key_arn = self.kms_backend.generate_data_key(
            key_id=key_id,
            encryption_context=encryption_context,
            number_of_bytes=number_of_bytes,
            key_spec=key_spec,
        )

        plaintext_response = base64.b64encode(plaintext).decode("utf-8")
        ciphertext_blob_response = base64.b64encode(ciphertext_blob).decode("utf-8")

        return json.dumps(
            {
                "CiphertextBlob": ciphertext_blob_response,
                "Plaintext": plaintext_response,
                "KeyId": key_arn,  # not alias
            }
        )

    def generate_data_key_without_plaintext(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_GenerateDataKeyWithoutPlaintext.html"""
        result = json.loads(self.generate_data_key())
        del result["Plaintext"]

        return json.dumps(result)

    def generate_random(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_GenerateRandom.html"""
        number_of_bytes = self._get_param("NumberOfBytes")

        if number_of_bytes and (number_of_bytes > 1024 or number_of_bytes < 1):
            raise ValidationException(
                (
                    "1 validation error detected: Value '{number_of_bytes:d}' at 'numberOfBytes' failed "
                    "to satisfy constraint: Member must have value less than or "
                    "equal to 1024"
                ).format(number_of_bytes=number_of_bytes)
            )

        entropy = os.urandom(number_of_bytes)

        response_entropy = base64.b64encode(entropy).decode("utf-8")

        return json.dumps({"Plaintext": response_entropy})

    def sign(self) -> str:
        key_id = self._get_param("KeyId")
        message = self._get_param("Message")
        message_type = self._get_param("MessageType")
        signing_algorithm = self._get_param("SigningAlgorithm")

        self._validate_key_id(key_id)

        if isinstance(message, str):
            message = message.encode("utf-8")

        if message == b"":
            raise ValidationException(
                "1 validation error detected: Value at 'Message' failed to satisfy constraint: Member must have length greater than or equal to 1"
            )

        if not message_type:
            message_type = "RAW"

        key_id, signature, signing_algorithm = self.kms_backend.sign(
            key_id=key_id,
            message=message,
            signing_algorithm=signing_algorithm,
        )

        signature_blob_response = base64.b64encode(signature).decode("utf-8")

        return json.dumps(
            {
                "KeyId": key_id,
                "Signature": signature_blob_response,
                "SigningAlgorithm": signing_algorithm,
            }
        )

    def verify(self) -> str:
        """https://docs.aws.amazon.com/kms/latest/APIReference/API_Verify.html"""
        key_id = self._get_param("KeyId")
        message = self._get_param("Message")
        message_type = self._get_param("MessageType")
        signature = self._get_param("Signature")
        signing_algorithm = self._get_param("SigningAlgorithm")

        self._validate_key_id(key_id)

        if not message_type:
            message_type = "RAW"

        if isinstance(message, str):
            message = message.encode("utf-8")

        if message == b"":
            raise ValidationException(
                "1 validation error detected: Value at 'Message' failed to satisfy constraint: Member must have length greater than or equal to 1"
            )

        if isinstance(signature, str):
            # we return base64 signatures, when signing
            signature = base64.b64decode(signature.encode("utf-8"))

        if signature == b"":
            raise ValidationException(
                "1 validation error detected: Value at 'Signature' failed to satisfy constraint: Member must have length greater than or equal to 1"
            )

        key_arn, signature_valid, signing_algorithm = self.kms_backend.verify(
            key_id=key_id,
            message=message,
            signature=signature,
            signing_algorithm=signing_algorithm,
        )

        return json.dumps(
            {
                "KeyId": key_arn,
                "SignatureValid": signature_valid,
                "SigningAlgorithm": signing_algorithm,
            }
        )

    def get_public_key(self) -> str:
        key_id = self._get_param("KeyId")

        self._validate_key_id(key_id)
        self._validate_cmk_id(key_id)
        key, public_key = self.kms_backend.get_public_key(key_id)
        return json.dumps(
            {
                "CustomerMasterKeySpec": key.key_spec,
                "EncryptionAlgorithms": key.encryption_algorithms,
                "KeyId": key.id,
                "KeyUsage": key.key_usage,
                "PublicKey": base64.b64encode(public_key).decode("UTF-8"),
                "SigningAlgorithms": key.signing_algorithms,
            }
        )


def _assert_default_policy(policy_name: str) -> None:
    if policy_name != "default":
        raise NotFoundException("No such policy exists")
