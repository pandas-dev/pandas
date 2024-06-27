import json
from typing import Any, Dict, List

from moto.core.responses import BaseResponse
from moto.secretsmanager.exceptions import (
    InvalidParameterException,
    InvalidRequestException,
    ValidationException,
)

from .models import SecretsManagerBackend, filter_keys, secretsmanager_backends


def _validate_filters(filters: List[Dict[str, Any]]) -> None:
    for idx, f in enumerate(filters):
        filter_key = f.get("Key", None)
        filter_values = f.get("Values", None)
        if filter_key is None:
            raise InvalidParameterException("Invalid filter key")
        if filter_key not in filter_keys():
            raise ValidationException(
                f"1 validation error detected: Value '{filter_key}' at 'filters.{(idx + 1)}.member.key' failed to satisfy constraint: "
                "Member must satisfy enum value set: [all, name, tag-key, description, tag-value]"
            )
        if filter_values is None:
            raise InvalidParameterException(
                f"Invalid filter values for key: {filter_key}"
            )


class SecretsManagerResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="secretsmanager")

    @property
    def backend(self) -> SecretsManagerBackend:
        return secretsmanager_backends[self.current_account][self.region]

    def cancel_rotate_secret(self) -> str:
        secret_id = self._get_param("SecretId")
        return self.backend.cancel_rotate_secret(secret_id=secret_id)

    def get_secret_value(self) -> str:
        secret_id = self._get_param("SecretId")
        version_id = self._get_param("VersionId")
        version_stage = self._get_param("VersionStage")
        value = self.backend.get_secret_value(
            secret_id=secret_id, version_id=version_id, version_stage=version_stage
        )
        return json.dumps(value)

    def batch_get_secret_value(self) -> str:
        secret_id_list = self._get_param("SecretIdList", if_none=[])
        filters = self._get_param("Filters", if_none=[])
        max_results = self._get_int_param("MaxResults")
        next_token = self._get_param("NextToken")

        secret_values, errors, next_token = self.backend.batch_get_secret_value(
            secret_id_list=secret_id_list,
            filters=filters,
            max_results=max_results,
            next_token=next_token,
        )
        return json.dumps(
            dict(SecretValues=secret_values, Errors=errors, NextToken=next_token)
        )

    def create_secret(self) -> str:
        name = self._get_param("Name")
        secret_string = self._get_param("SecretString")
        secret_binary = self._get_param("SecretBinary")
        description = self._get_param("Description", if_none="")
        tags = self._get_param("Tags", if_none=[])
        kms_key_id = self._get_param("KmsKeyId")
        client_request_token = self._get_param("ClientRequestToken")
        replica_regions = self._get_param("AddReplicaRegions", [])
        force_overwrite = self._get_bool_param("ForceOverwriteReplicaSecret", False)

        return self.backend.create_secret(
            name=name,
            secret_string=secret_string,
            secret_binary=secret_binary,
            description=description,
            tags=tags,
            kms_key_id=kms_key_id,
            client_request_token=client_request_token,
            replica_regions=replica_regions,
            force_overwrite=force_overwrite,
        )

    def update_secret(self) -> str:
        secret_id = self._get_param("SecretId")
        secret_string = self._get_param("SecretString")
        secret_binary = self._get_param("SecretBinary")
        description = self._get_param("Description")
        client_request_token = self._get_param("ClientRequestToken")
        kms_key_id = self._get_param("KmsKeyId", if_none=None)
        return self.backend.update_secret(
            secret_id=secret_id,
            secret_string=secret_string,
            secret_binary=secret_binary,
            client_request_token=client_request_token,
            kms_key_id=kms_key_id,
            description=description,
        )

    def get_random_password(self) -> str:
        password_length = self._get_param("PasswordLength", if_none=32)
        exclude_characters = self._get_param("ExcludeCharacters", if_none="")
        exclude_numbers = self._get_param("ExcludeNumbers", if_none=False)
        exclude_punctuation = self._get_param("ExcludePunctuation", if_none=False)
        exclude_uppercase = self._get_param("ExcludeUppercase", if_none=False)
        exclude_lowercase = self._get_param("ExcludeLowercase", if_none=False)
        include_space = self._get_param("IncludeSpace", if_none=False)
        require_each_included_type = self._get_param(
            "RequireEachIncludedType", if_none=True
        )
        return self.backend.get_random_password(
            password_length=password_length,
            exclude_characters=exclude_characters,
            exclude_numbers=exclude_numbers,
            exclude_punctuation=exclude_punctuation,
            exclude_uppercase=exclude_uppercase,
            exclude_lowercase=exclude_lowercase,
            include_space=include_space,
            require_each_included_type=require_each_included_type,
        )

    def describe_secret(self) -> str:
        secret_id = self._get_param("SecretId")
        secret = self.backend.describe_secret(secret_id=secret_id)
        return json.dumps(secret.to_dict())

    def rotate_secret(self) -> str:
        client_request_token = self._get_param("ClientRequestToken")
        rotation_lambda_arn = self._get_param("RotationLambdaARN")
        rotation_rules = self._get_param("RotationRules")
        secret_id = self._get_param("SecretId")
        rotate_immediately = self._get_bool_param("RotateImmediately", True)
        return self.backend.rotate_secret(
            secret_id=secret_id,
            client_request_token=client_request_token,
            rotation_lambda_arn=rotation_lambda_arn,
            rotation_rules=rotation_rules,
            rotate_immediately=rotate_immediately,
        )

    def put_secret_value(self) -> str:
        secret_id = self._get_param("SecretId", if_none="")
        secret_string = self._get_param("SecretString")
        secret_binary = self._get_param("SecretBinary")
        client_request_token = self._get_param("ClientRequestToken")
        if not secret_binary and not secret_string:
            raise InvalidRequestException(
                "You must provide either SecretString or SecretBinary."
            )
        version_stages = self._get_param("VersionStages", if_none=["AWSCURRENT"])
        if not isinstance(version_stages, list):
            version_stages = [version_stages]

        return self.backend.put_secret_value(
            secret_id=secret_id,
            secret_binary=secret_binary,
            secret_string=secret_string,
            version_stages=version_stages,
            client_request_token=client_request_token,
        )

    def list_secret_version_ids(self) -> str:
        secret_id = self._get_param("SecretId", if_none="")
        return self.backend.list_secret_version_ids(secret_id=secret_id)

    def list_secrets(self) -> str:
        filters = self._get_param("Filters", if_none=[])
        _validate_filters(filters)
        max_results = self._get_int_param("MaxResults")
        next_token = self._get_param("NextToken")
        secret_list, next_token = self.backend.list_secrets(
            filters=filters, max_results=max_results, next_token=next_token
        )
        return json.dumps(dict(SecretList=secret_list, NextToken=next_token))

    def delete_secret(self) -> str:
        secret_id = self._get_param("SecretId")
        recovery_window_in_days = self._get_param("RecoveryWindowInDays")
        force_delete_without_recovery = self._get_param("ForceDeleteWithoutRecovery")
        arn, name, deletion_date = self.backend.delete_secret(
            secret_id=secret_id,
            recovery_window_in_days=recovery_window_in_days,
            force_delete_without_recovery=force_delete_without_recovery,
        )
        return json.dumps(dict(ARN=arn, Name=name, DeletionDate=deletion_date))

    def restore_secret(self) -> str:
        secret_id = self._get_param("SecretId")
        arn, name = self.backend.restore_secret(secret_id=secret_id)
        return json.dumps(dict(ARN=arn, Name=name))

    def get_resource_policy(self) -> str:
        secret_id = self._get_param("SecretId")
        return self.backend.get_resource_policy(secret_id=secret_id)

    def put_resource_policy(self) -> str:
        secret_id = self._get_param("SecretId")
        policy = self._get_param("ResourcePolicy")
        arn, name = self.backend.put_resource_policy(secret_id, policy)
        return json.dumps(dict(ARN=arn, Name=name))

    def delete_resource_policy(self) -> str:
        secret_id = self._get_param("SecretId")
        arn, name = self.backend.delete_resource_policy(secret_id)
        return json.dumps(dict(ARN=arn, Name=name))

    def tag_resource(self) -> str:
        secret_id = self._get_param("SecretId")
        tags = self._get_param("Tags", if_none=[])
        self.backend.tag_resource(secret_id, tags)
        return "{}"

    def untag_resource(self) -> str:
        secret_id = self._get_param("SecretId")
        tag_keys = self._get_param("TagKeys", if_none=[])
        self.backend.untag_resource(secret_id=secret_id, tag_keys=tag_keys)
        return "{}"

    def update_secret_version_stage(self) -> str:
        secret_id = self._get_param("SecretId")
        version_stage = self._get_param("VersionStage")
        remove_from_version_id = self._get_param("RemoveFromVersionId")
        move_to_version_id = self._get_param("MoveToVersionId")
        arn, name = self.backend.update_secret_version_stage(
            secret_id=secret_id,
            version_stage=version_stage,
            remove_from_version_id=remove_from_version_id,
            move_to_version_id=move_to_version_id,
        )
        return json.dumps({"ARN": arn, "Name": name})

    def replicate_secret_to_regions(self) -> str:
        secret_id = self._get_param("SecretId")
        replica_regions = self._get_param("AddReplicaRegions", [])
        force_overwrite = self._get_bool_param("ForceOverwriteReplicaSecret", False)

        arn, statuses = self.backend.replicate_secret_to_regions(
            secret_id, replica_regions, force_overwrite=force_overwrite
        )
        return json.dumps({"ARN": arn, "ReplicationStatus": statuses})

    def remove_regions_from_replication(self) -> str:
        secret_id = self._get_param("SecretId")
        replica_regions = self._get_param("RemoveReplicaRegions", [])

        arn, statuses = self.backend.remove_regions_from_replication(
            secret_id, replica_regions
        )
        return json.dumps({"ARN": arn, "ReplicationStatus": statuses})
