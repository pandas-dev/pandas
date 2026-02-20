import datetime
import json
import time
from typing import Any, Optional, Union

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcfromtimestamp, utcnow
from moto.moto_api._internal import mock_random

from .exceptions import (
    ClientError,
    InvalidParameterException,
    InvalidRequestException,
    OperationNotPermittedOnReplica,
    ResourceExistsException,
    ResourceNotFoundException,
    SecretHasNoValueException,
    SecretNotFoundException,
    SecretStageVersionMismatchException,
)
from .list_secrets.filters import (
    description_filter,
    filter_all,
    name_filter,
    owning_service_filter,
    tag_key,
    tag_value,
)
from .utils import (
    SecretsManagerSecretIdentifier,
    get_secret_name_from_partial_arn,
    random_password,
)

MAX_RESULTS_DEFAULT = 100


def filter_primary_region(
    secret: Union["FakeSecret", "ReplicaSecret"], values: list[str]
) -> bool:
    if isinstance(secret, FakeSecret):
        return len(secret.replicas) > 0 and secret.region in values
    elif isinstance(secret, ReplicaSecret):
        return secret.source.region in values


_filter_functions = {
    "all": filter_all,
    "name": name_filter,
    "description": description_filter,
    "tag-key": tag_key,
    "tag-value": tag_value,
    "primary-region": filter_primary_region,
    "owning-service": owning_service_filter,
}


def filter_keys() -> list[str]:
    return list(_filter_functions.keys())


def _matches(
    secret: Union["FakeSecret", "ReplicaSecret"], filters: list[dict[str, Any]]
) -> bool:
    is_match = True

    for f in filters:
        # Filter names are pre-validated in the resource layer
        filter_function = _filter_functions.get(f["Key"])
        is_match = is_match and filter_function(secret, f["Values"])  # type: ignore

    return is_match


class FakeSecret(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        secret_id: str,
        secret_version: dict[str, Any],
        version_id: str,
        secret_string: Optional[str] = None,
        secret_binary: Optional[str] = None,
        description: Optional[str] = None,
        tags: Optional[list[dict[str, str]]] = None,
        kms_key_id: Optional[str] = None,
        version_stages: Optional[list[str]] = None,
        last_changed_date: Optional[int] = None,
        created_date: Optional[int] = None,
        replica_regions: Optional[list[dict[str, str]]] = None,
        force_overwrite: bool = False,
    ):
        self.secret_id = secret_id
        self.name = secret_id
        self.arn = SecretsManagerSecretIdentifier(
            account_id, region_name, secret_id
        ).generate(tags=tags)
        self.account_id = account_id
        self.region = region_name
        self.secret_string = secret_string
        self.secret_binary = secret_binary
        self.description = description
        self.tags = tags or None
        self.kms_key_id = kms_key_id
        self.version_stages = version_stages
        self.last_changed_date = last_changed_date
        self.created_date = created_date
        # We should only return Rotation details after it's been requested
        self.rotation_requested = False
        self.rotation_enabled = False
        self.rotation_lambda_arn = ""
        self.auto_rotate_after_days = 0
        self.deleted_date: Optional[float] = None
        self.policy: Optional[str] = None
        self.next_rotation_date: Optional[int] = None
        self.last_rotation_date: Optional[int] = None

        self.versions: dict[str, dict[str, Any]] = {}
        if secret_string or secret_binary:
            self.versions = {version_id: secret_version}
            self.set_default_version_id(version_id)
        else:
            self.set_default_version_id(None)

        self.replicas = self.create_replicas(
            replica_regions or [], force_overwrite=force_overwrite
        )

    @property
    def owning_service(self) -> Optional[str]:
        for tag in self.tags or []:
            if tag["Key"] == "aws:secretsmanager:owningService":
                return tag["Value"]
        return None

    def create_replicas(
        self, replica_regions: list[dict[str, str]], force_overwrite: bool
    ) -> list["ReplicaSecret"]:
        # Validate first, before we create anything
        for replica_config in replica_regions or []:
            if replica_config["Region"] == self.region:
                raise InvalidParameterException("Invalid replica region.")

        replicas: list[ReplicaSecret] = []
        for replica_config in replica_regions or []:
            replica_region = replica_config["Region"]
            backend = secretsmanager_backends[self.account_id][replica_region]
            if self.name in backend.secrets:
                if force_overwrite:
                    backend.secrets.pop(self.name)
                    replica = ReplicaSecret(self, replica_region)
                    backend.secrets[replica.arn] = replica
                else:
                    message = f"Replication failed: Secret name simple already exists in region {backend.region_name}."
                    replica = ReplicaSecret(self, replica_region, "Failed", message)
            else:
                replica = ReplicaSecret(self, replica_region)
                backend.secrets[replica.arn] = replica
            replicas.append(replica)
        return replicas

    def update(
        self,
        description: Optional[str] = None,
        tags: Optional[list[dict[str, str]]] = None,
        kms_key_id: Optional[str] = None,
        last_changed_date: Optional[int] = None,
    ) -> None:
        self.description = description
        self.tags = tags or None
        if last_changed_date is not None:
            self.last_changed_date = last_changed_date

        if kms_key_id is not None:
            self.kms_key_id = kms_key_id

    def set_default_version_id(self, version_id: Optional[str]) -> None:
        self.default_version_id = version_id

    def reset_default_version(
        self, secret_version: dict[str, Any], version_id: str
    ) -> None:
        # remove all old AWSPREVIOUS stages
        for old_version in self.versions.values():
            if "AWSPREVIOUS" in old_version["version_stages"]:
                old_version["version_stages"].remove("AWSPREVIOUS")

        if self.default_version_id:
            # set old AWSCURRENT secret to AWSPREVIOUS
            previous_current_version_id = self.default_version_id
            self.versions[previous_current_version_id]["version_stages"] = [
                "AWSPREVIOUS"
            ]

        self.versions[version_id] = secret_version
        self.default_version_id = version_id

    def remove_version_stages_from_old_versions(
        self, version_stages: list[str]
    ) -> None:
        for version_stage in version_stages:
            for old_version in self.versions.values():
                if version_stage in old_version["version_stages"]:
                    old_version["version_stages"].remove(version_stage)

    def delete(self, deleted_date: float) -> None:
        self.deleted_date = deleted_date

    def restore(self) -> None:
        self.deleted_date = None

    def is_deleted(self) -> bool:
        return self.deleted_date is not None

    def to_short_dict(
        self,
        include_version_stages: bool = False,
        version_id: Optional[str] = None,
        include_version_id: bool = True,
    ) -> str:
        if not version_id:
            version_id = self.default_version_id
        dct: dict[str, Any] = {
            "ARN": self.arn,
            "Name": self.name,
        }
        if include_version_id and version_id:
            dct["VersionId"] = version_id
        if version_id and include_version_stages:
            dct["VersionStages"] = self.versions[version_id]["version_stages"]
        if self.replicas:
            dct["ReplicationStatus"] = [replica.config for replica in self.replicas]
        return json.dumps(dct)

    def to_dict(self) -> dict[str, Any]:
        version_id_to_stages = self._form_version_ids_to_stages()

        dct: dict[str, Any] = {
            "ARN": self.arn,
            "Name": self.name,
            "LastChangedDate": self.last_changed_date,
            "LastAccessedDate": None,
            "NextRotationDate": self.next_rotation_date,
            "DeletedDate": self.deleted_date,
            "CreatedDate": self.created_date,
        }
        if self.kms_key_id != SecretsManagerBackend.DEFAULT_KMS_KEY_ALIAS:
            dct["KmsKeyId"] = self.kms_key_id
        if self.owning_service is not None:
            dct["OwningService"] = self.owning_service
        if self.tags is not None:
            dct["Tags"] = self.tags
        if self.description:
            dct["Description"] = self.description
        if self.versions:
            dct.update(
                {
                    # Key used by describe_secret
                    "VersionIdsToStages": version_id_to_stages,
                    # Key used by list_secrets
                    "SecretVersionsToStages": version_id_to_stages,
                }
            )
        if self.rotation_requested:
            dct.update(
                {
                    "RotationEnabled": self.rotation_enabled,
                    "RotationLambdaARN": self.rotation_lambda_arn,
                    "RotationRules": {
                        "AutomaticallyAfterDays": self.auto_rotate_after_days
                    },
                    "LastRotatedDate": self.last_rotation_date,
                }
            )
        if self.replicas:
            dct["ReplicationStatus"] = [replica.config for replica in self.replicas]
        return dct

    def _form_version_ids_to_stages(self) -> dict[str, str]:
        version_id_to_stages = {}
        for key, value in self.versions.items():
            version_id_to_stages[key] = value["version_stages"]

        return version_id_to_stages


class ReplicaSecret:
    def __init__(
        self,
        source: FakeSecret,
        region: str,
        status: Optional[str] = None,
        message: Optional[str] = None,
    ):
        self.source = source
        self.arn = source.arn.replace(source.region, region)
        self.region = region
        self.status = status or "InSync"
        self.message = message or "Replication succeeded"
        self.has_replica = status is None
        self.config = {
            "Region": self.region,
            "KmsKeyId": SecretsManagerBackend.DEFAULT_KMS_KEY_ALIAS,
            "Status": self.status,
            "StatusMessage": self.message,
        }

    def is_deleted(self) -> bool:
        return False

    def to_dict(self) -> dict[str, Any]:
        dct = self.source.to_dict()
        dct["ARN"] = self.arn
        dct["PrimaryRegion"] = self.source.region
        return dct

    @property
    def default_version_id(self) -> Optional[str]:
        return self.source.default_version_id

    @property
    def versions(self) -> dict[str, dict[str, Any]]:  # type: ignore[misc]
        return self.source.versions

    @property
    def name(self) -> str:
        return self.source.name

    @property
    def secret_id(self) -> str:
        return self.source.secret_id

    @property
    def description(self) -> Optional[str]:
        return self.source.description

    @property
    def tags(self) -> Optional[list[dict[str, str]]]:
        return self.source.tags

    @property
    def owning_service(self) -> Optional[str]:
        return self.source.owning_service


class SecretsStore(dict[str, Union[FakeSecret, ReplicaSecret]]):
    # Parameters to this dictionary can be three possible values:
    # names, full ARNs, and partial ARNs
    # Every retrieval method should check which type of input it receives

    def __setitem__(self, key: str, value: Union[FakeSecret, ReplicaSecret]) -> None:
        super().__setitem__(key, value)

    def __getitem__(self, key: str) -> Union[FakeSecret, ReplicaSecret]:
        for secret in dict.values(self):
            if secret.arn == key or secret.name == key:
                return secret
        name = get_secret_name_from_partial_arn(key)
        return super().__getitem__(name)

    def __contains__(self, key: str) -> bool:  # type: ignore
        for secret in dict.values(self):
            if secret.arn == key or secret.name == key:
                return True
        name = get_secret_name_from_partial_arn(key)
        return dict.__contains__(self, name)  # type: ignore

    def get(self, key: str) -> Optional[Union[FakeSecret, ReplicaSecret]]:  # type: ignore
        for secret in dict.values(self):
            if secret.arn == key or secret.name == key:
                return secret
        name = get_secret_name_from_partial_arn(key)
        return super().get(name)

    def pop(self, key: str) -> Optional[Union[FakeSecret, ReplicaSecret]]:  # type: ignore
        for secret in dict.values(self):
            if secret.arn == key or secret.name == key:
                key = secret.name
        name = get_secret_name_from_partial_arn(key)
        return super().pop(name, None)


class SecretsManagerBackend(BaseBackend):
    DEFAULT_KMS_KEY_ALIAS = "alias/aws/secretsmanager"

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.secrets = SecretsStore()

    def _is_valid_identifier(self, identifier: str) -> bool:
        return identifier in self.secrets

    def _unix_time_secs(self, dt: datetime.datetime) -> float:
        epoch = utcfromtimestamp(0)
        return (dt - epoch).total_seconds()

    def _client_request_token_validator(self, client_request_token: str) -> None:
        token_length = len(client_request_token)
        if token_length < 32 or token_length > 64:
            msg = "ClientRequestToken must be 32-64 characters long."
            raise InvalidParameterException(msg)

    def _from_client_request_token(self, client_request_token: Optional[str]) -> str:
        if client_request_token:
            self._client_request_token_validator(client_request_token)
            return client_request_token
        else:
            return str(mock_random.uuid4())

    def cancel_rotate_secret(self, secret_id: str) -> str:
        if not self._is_valid_identifier(secret_id):
            raise SecretNotFoundException()

        secret = self.secrets[secret_id]
        if isinstance(secret, ReplicaSecret):
            raise OperationNotPermittedOnReplica
        if secret.is_deleted():
            raise InvalidRequestException(
                "You tried to perform the operation on a secret that's currently marked deleted."
            )

        if not secret.rotation_lambda_arn:
            # This response doesn't make much sense for  `CancelRotateSecret`, but this is what AWS has documented ...
            # https://docs.aws.amazon.com/secretsmanager/latest/apireference/API_CancelRotateSecret.html
            raise InvalidRequestException(
                "You tried to enable rotation on a secret that doesn't already have a Lambda function ARN configured"
                "and you didn't include such an ARN as a parameter in this call."
            )

        secret.rotation_enabled = False
        return secret.to_short_dict()

    def get_secret_value(
        self, secret_id: str, version_id: str, version_stage: str
    ) -> dict[str, Any]:
        if not self._is_valid_identifier(secret_id):
            raise SecretNotFoundException()

        if version_id and version_stage:
            versions_dict = self.secrets[secret_id].versions
            if (
                version_id in versions_dict
                and version_stage not in versions_dict[version_id]["version_stages"]
            ):
                raise SecretStageVersionMismatchException()

        version_id_provided = version_id is not None
        if not version_id and version_stage:
            # set version_id to match version_stage
            versions_dict = self.secrets[secret_id].versions
            for ver_id, ver_val in versions_dict.items():
                if version_stage in ver_val["version_stages"]:
                    version_id = ver_id
                    break
            if not version_id:
                raise SecretNotFoundException()

        # TODO check this part
        if self.secrets[secret_id].is_deleted():
            raise InvalidRequestException(
                "An error occurred (InvalidRequestException) when calling the GetSecretValue operation: You tried to \
                perform the operation on a secret that's currently marked deleted."
            )

        secret = self.secrets[secret_id]
        version_id = version_id or secret.default_version_id or "AWSCURRENT"

        secret_version = secret.versions.get(version_id)
        if not secret_version:
            _type = "staging label" if not version_id_provided else "VersionId"
            raise ResourceNotFoundException(
                f"Secrets Manager can't find the specified secret value for {_type}: {version_id}"
            )

        response_data = {
            "ARN": secret.arn,
            "Name": secret.name,
            "VersionId": secret_version["version_id"],
            "VersionStages": secret_version["version_stages"],
            "CreatedDate": secret_version["createdate"],
        }

        if "secret_string" in secret_version:
            response_data["SecretString"] = secret_version["secret_string"]

        if "secret_binary" in secret_version:
            response_data["SecretBinary"] = secret_version["secret_binary"]

        if (
            "secret_string" not in secret_version
            and "secret_binary" not in secret_version
        ):
            raise SecretHasNoValueException(version_stage or "AWSCURRENT")

        return response_data

    def batch_get_secret_value(
        self,
        secret_id_list: Optional[list[str]] = None,
        filters: Optional[list[dict[str, Any]]] = None,
        max_results: Optional[int] = None,
        next_token: Optional[str] = None,
    ) -> tuple[list[dict[str, Any]], list[Any], Optional[str]]:
        secret_list = []
        errors: list[Any] = []
        if secret_id_list and filters:
            raise InvalidParameterException(
                "Either 'SecretIdList' or 'Filters' must be provided, but not both."
            )

        if max_results and not filters:
            raise InvalidParameterException(
                "'Filters' not specified. 'Filters' must also be specified when 'MaxResults' is provided."
            )

        if secret_id_list:
            for secret_id in secret_id_list:
                # TODO perhaps there should be a check if the secret id is valid identifier
                # and add an error to the list if not
                try:
                    secret_list.append(self.get_secret_value(secret_id, "", ""))
                except (SecretNotFoundException, InvalidRequestException) as e:
                    errors.append(
                        {
                            "SecretId": secret_id,
                            "ErrorCode": e.error_type,
                            "Message": e.message,
                        }
                    )

        if filters:
            for secret in self.secrets.values():
                if _matches(secret, filters):
                    if isinstance(secret, FakeSecret):
                        secret_list.append(
                            self.get_secret_value(secret.secret_id, "", "")
                        )
                    elif isinstance(secret, ReplicaSecret):
                        secret_list.append(
                            self.get_secret_value(secret.source.secret_id, "", "")
                        )

        secret_page, new_next_token = self._get_secret_values_page_and_next_token(
            secret_list, max_results, next_token
        )
        return secret_page, errors, new_next_token

    def update_secret(
        self,
        secret_id: str,
        secret_string: Optional[str] = None,
        secret_binary: Optional[str] = None,
        client_request_token: Optional[str] = None,
        kms_key_id: Optional[str] = None,
        description: Optional[str] = None,
    ) -> str:
        # error if secret does not exist
        if secret_id not in self.secrets:
            raise SecretNotFoundException()

        existing_secret_info = self._check_with_existing_secrets_and_versions(
            secret_id, client_request_token, secret_string, secret_binary
        )
        if existing_secret_info:
            return existing_secret_info

        secret = self.secrets[secret_id]
        if secret.is_deleted():
            raise InvalidRequestException(
                "An error occurred (InvalidRequestException) when calling the UpdateSecret operation: "
                "You can't perform this operation on the secret because it was marked for deletion."
            )

        tags = secret.tags
        description = description or secret.description

        secret, new_version = self._add_secret(
            secret_id,
            secret_string=secret_string,
            secret_binary=secret_binary,
            description=description,
            version_id=client_request_token,
            tags=tags,
            kms_key_id=kms_key_id,
        )

        return secret.to_short_dict(include_version_id=new_version)

    def create_secret(
        self,
        name: str,
        secret_string: Optional[str],
        secret_binary: Optional[str],
        description: Optional[str],
        tags: Optional[list[dict[str, str]]],
        kms_key_id: Optional[str],
        client_request_token: Optional[str],
        replica_regions: list[dict[str, str]],
        force_overwrite: bool,
    ) -> str:
        existing_secret = self._check_with_existing_secrets_and_versions(
            name, client_request_token, secret_string, secret_binary
        )
        if existing_secret:
            return existing_secret

        secret, new_version = self._add_secret(
            name,
            secret_string=secret_string,
            secret_binary=secret_binary,
            description=description,
            tags=tags,
            kms_key_id=kms_key_id,
            version_id=client_request_token,
            replica_regions=replica_regions,
            force_overwrite=force_overwrite,
        )

        return secret.to_short_dict(include_version_id=new_version)

    def _check_with_existing_secrets_and_versions(
        self,
        secret_name: str,
        client_request_token: Optional[str],
        secret_string: Optional[str],
        secret_binary: Optional[str],
    ) -> Optional[str]:
        """
        Check if a secret with the given name and version ID already exists.
        If they do and the original operation is intending to modify it, that will be flagged.
        Since we can only add new versions of secrets, not modify existing ones
        """
        if secret_name not in self.secrets.keys():
            # Nothing to validate here
            return None
        existing_secret = self.secrets[secret_name]
        if isinstance(existing_secret, ReplicaSecret):
            raise OperationNotPermittedOnReplica
        # If the secret already exists, and a client request token is provided,
        # we need to check if it matches the existing secret's version ID.
        if not client_request_token:
            # No version identifier was provided to compare with
            raise ResourceExistsException(
                "A resource with the ID you requested already exists."
            )
        # Check if client_request_token is part of any version of this secret
        if client_request_token not in existing_secret.versions:
            return None
        # Check if the secret_string/secret_binary values corresponding to this version matches that of the current request
        matching_secret_version = existing_secret.versions[client_request_token]
        # If this version stage label is AWSPENDING, it means this version was created as part of rotation
        if "AWSPENDING" in matching_secret_version.get("version_stages", []):
            return None
        if (
            matching_secret_version.get("secret_string") == secret_string
            and matching_secret_version.get("secret_binary") == secret_binary
        ):
            # If they match, we can return the existing secret without error
            return existing_secret.to_short_dict(
                include_version_id=True, version_id=client_request_token
            )
        # If they do not match, though, then the request fails since we cannot modify
        # an existing version
        raise ResourceExistsException(
            f"You can't use ClientRequestToken {client_request_token} because that value is already in use for a version of secret {existing_secret.arn}"
        )

    def _add_secret(
        self,
        secret_id: str,
        secret_string: Optional[str] = None,
        secret_binary: Optional[str] = None,
        description: Optional[str] = None,
        tags: Optional[list[dict[str, str]]] = None,
        kms_key_id: Optional[str] = None,
        version_id: Optional[str] = None,
        version_stages: Optional[list[str]] = None,
        replica_regions: Optional[list[dict[str, str]]] = None,
        force_overwrite: bool = False,
    ) -> tuple[FakeSecret, bool]:
        if version_stages is None:
            version_stages = ["AWSCURRENT"]

        version_id = self._from_client_request_token(version_id)

        secret_version = {
            "createdate": int(time.time()),
            "version_id": version_id,
            "version_stages": version_stages,
        }
        if secret_string is not None:
            secret_version["secret_string"] = secret_string

        if secret_binary is not None:
            secret_version["secret_binary"] = secret_binary

        new_version = secret_string is not None or secret_binary is not None

        update_time = int(time.time())
        if secret_id in self.secrets:
            secret = self.secrets[secret_id]
            if isinstance(secret, ReplicaSecret):
                raise OperationNotPermittedOnReplica

            secret.update(description, tags, kms_key_id, last_changed_date=update_time)

            if new_version:
                if "AWSCURRENT" in version_stages:
                    secret.reset_default_version(secret_version, version_id)
                else:
                    secret.remove_version_stages_from_old_versions(version_stages)
                    secret.versions[version_id] = secret_version
        else:
            secret = FakeSecret(
                account_id=self.account_id,
                region_name=self.region_name,
                secret_id=secret_id,
                secret_string=secret_string,
                secret_binary=secret_binary,
                description=description,
                tags=tags,
                kms_key_id=kms_key_id,
                last_changed_date=update_time,
                created_date=update_time,
                version_id=version_id,
                secret_version=secret_version,
                replica_regions=replica_regions,
                force_overwrite=force_overwrite,
            )
            self.secrets[secret_id] = secret

        return secret, new_version

    def create_managed_secret(
        self,
        service_name: str,
        secret_id: str,
        secret_string: Optional[str] = None,
        secret_binary: Optional[str] = None,
        description: Optional[str] = None,
        tags: Optional[list[dict[str, str]]] = None,
        kms_key_id: Optional[str] = None,
        version_id: Optional[str] = None,
        replica_regions: Optional[list[dict[str, str]]] = None,
        force_overwrite: bool = False,
    ) -> FakeSecret:
        """Create an AWS managed secret for the specified service name."""
        if kms_key_id is None:
            kms_key_id = self.DEFAULT_KMS_KEY_ALIAS
        managed_tag = {
            "Key": "aws:secretsmanager:owningService",
            "Value": service_name,
        }
        if tags is None:
            tags = [managed_tag]
        else:
            tags.append(managed_tag)
        secret, _ = self._add_secret(
            secret_id,
            secret_string=secret_string,
            secret_binary=secret_binary,
            description=description,
            tags=tags,
            kms_key_id=kms_key_id,
            version_id=version_id,
            replica_regions=replica_regions,
            force_overwrite=force_overwrite,
        )
        return secret

    def put_secret_value(
        self,
        secret_id: str,
        secret_string: str,
        secret_binary: str,
        client_request_token: str,
        version_stages: list[str],
    ) -> str:
        if not self._is_valid_identifier(secret_id):
            raise SecretNotFoundException()
        else:
            secret = self.secrets[secret_id]
            if isinstance(secret, ReplicaSecret):
                raise OperationNotPermittedOnReplica
            tags = secret.tags
            description = secret.description

        version_id = self._from_client_request_token(client_request_token)
        existing_secret = self._check_with_existing_secrets_and_versions(
            secret.name, version_id, secret_string, secret_binary
        )
        # If it exists, then return the existing secret
        if existing_secret:
            return existing_secret

        # If it is the first version add AWSCURRENT to the versions
        if not secret.versions and "AWSCURRENT" not in version_stages:
            version_stages.append("AWSCURRENT")

        secret, _ = self._add_secret(
            secret_id,
            secret_string,
            secret_binary,
            version_id=version_id,
            description=description,
            tags=tags,
            version_stages=version_stages,
        )

        return secret.to_short_dict(include_version_stages=True, version_id=version_id)

    def describe_secret(self, secret_id: str) -> Union[FakeSecret, ReplicaSecret]:
        if not self._is_valid_identifier(secret_id):
            raise SecretNotFoundException()

        return self.secrets[secret_id]

    def rotate_secret(
        self,
        secret_id: str,
        client_request_token: Optional[str] = None,
        rotation_lambda_arn: Optional[str] = None,
        rotation_rules: Optional[dict[str, Any]] = None,
        rotate_immediately: bool = True,
    ) -> str:
        rotation_days = "AutomaticallyAfterDays"

        if not self._is_valid_identifier(secret_id):
            raise SecretNotFoundException()

        secret = self.secrets[secret_id]
        if isinstance(secret, ReplicaSecret):
            raise OperationNotPermittedOnReplica
        if secret.is_deleted():
            raise InvalidRequestException(
                "An error occurred (InvalidRequestException) when calling the RotateSecret operation: You tried to \
                perform the operation on a secret that's currently marked deleted."
            )

        if rotation_lambda_arn:
            if len(rotation_lambda_arn) > 2048:
                msg = "RotationLambdaARN must <= 2048 characters long."
                raise InvalidParameterException(msg)

        if rotation_rules:
            if rotation_days in rotation_rules:
                rotation_period = rotation_rules[rotation_days]
                if rotation_period < 1 or rotation_period > 1000:
                    msg = "RotationRules.AutomaticallyAfterDays must be within 1-1000."
                    raise InvalidParameterException(msg)

                secret.next_rotation_date = int(time.time()) + (
                    int(rotation_period) * 86400
                )

        # The rotation function must end with the versions of the secret in
        # one of two states:
        #
        #  - The AWSPENDING and AWSCURRENT staging labels are attached to the
        #    same version of the secret, or
        #  - The AWSPENDING staging label is not attached to any version of the secret.
        #
        # If the AWSPENDING staging label is present but not attached to the same
        # version as AWSCURRENT then any later invocation of RotateSecret assumes
        # that a previous rotation request is still in progress and returns an error.
        try:
            version = next(
                version
                for version in secret.versions.values()
                if "AWSPENDING" in version["version_stages"]
            )
            if "AWSCURRENT" in version["version_stages"]:
                msg = "Previous rotation request is still in progress."
                raise InvalidRequestException(msg)

        except StopIteration:
            # Pending is not present in any version
            pass

        if secret.versions:
            if client_request_token:
                self._client_request_token_validator(client_request_token)
                new_version_id = client_request_token
            else:
                new_version_id = str(mock_random.uuid4())

            # We add a "pending" stage. The previous version remains as "current" for now.
            # Caller is responsible for creating the new secret in the Lambda
            secret_version = {
                "createdate": int(time.time()),
                "version_id": new_version_id,
                "version_stages": ["AWSPENDING"],
            }
            if not rotate_immediately:
                if secret.secret_string is not None:
                    secret_version["secret_string"] = secret.secret_string
                if secret.secret_binary is not None:
                    secret_version["secret_binary"] = secret.secret_binary

            secret.remove_version_stages_from_old_versions(["AWSPENDING"])
            secret.versions[new_version_id] = secret_version

        secret.rotation_requested = True
        secret.rotation_lambda_arn = rotation_lambda_arn or ""
        if rotation_rules:
            secret.auto_rotate_after_days = rotation_rules.get(rotation_days, 0)
        if secret.auto_rotate_after_days > 0:
            secret.rotation_enabled = True

        # Begin the rotation process for the given secret by invoking the lambda function.
        if secret.rotation_lambda_arn:
            from moto.awslambda.utils import get_backend

            lambda_backend = get_backend(self.account_id, self.region_name)

            request_headers: dict[str, Any] = {}
            response_headers: dict[str, Any] = {}

            try:
                lambda_backend.get_function(secret.rotation_lambda_arn)
            except Exception:
                msg = f"Resource not found for ARN '{secret.rotation_lambda_arn}'."
                raise ResourceNotFoundException(msg)

            rotation_steps = ["create", "set", "test", "finish"]
            if not rotate_immediately:
                # if you don't immediately rotate the secret,
                # Secrets Manager tests the rotation configuration by running the testSecretstep of the Lambda rotation function.
                rotation_steps = ["test"]
            for step in rotation_steps:
                lambda_backend.invoke(
                    secret.rotation_lambda_arn,
                    qualifier=None,
                    body=json.dumps(
                        {
                            "Step": step + "Secret",
                            "SecretId": secret.name,
                            "ClientRequestToken": new_version_id,
                        }
                    ),
                    headers=request_headers,
                    response_headers=response_headers,
                )
            if rotate_immediately:
                # If we don't rotate, we only invoke the testSecret step
                # This should be done with the existing (old) version ID
                secret.set_default_version_id(new_version_id)

        elif secret.versions:
            # AWS will always require a Lambda ARN
            # without that, Moto can still apply the 'AWSCURRENT'-label
            # This only makes sense if we have a version
            secret.reset_default_version(
                secret.versions[new_version_id], new_version_id
            )
            secret.versions[new_version_id]["version_stages"] = ["AWSCURRENT"]

        secret.last_rotation_date = int(time.time())
        return secret.to_short_dict()

    def get_random_password(
        self,
        password_length: int,
        exclude_characters: str,
        exclude_numbers: bool,
        exclude_punctuation: bool,
        exclude_uppercase: bool,
        exclude_lowercase: bool,
        include_space: bool,
        require_each_included_type: bool,
    ) -> str:
        # password size must have value less than or equal to 4096
        if password_length > 4096:
            raise ClientError(
                f"ClientError: An error occurred (ValidationException) \
                when calling the GetRandomPassword operation: 1 validation error detected: Value '{password_length}' at 'passwordLength' \
                failed to satisfy constraint: Member must have value less than or equal to 4096"
            )
        if password_length < 4:
            raise InvalidParameterException(
                "InvalidParameterException: An error occurred (InvalidParameterException) \
                when calling the GetRandomPassword operation: Password length is too short based on the required types."
            )

        return json.dumps(
            {
                "RandomPassword": random_password(
                    password_length,
                    exclude_characters,
                    exclude_numbers,
                    exclude_punctuation,
                    exclude_uppercase,
                    exclude_lowercase,
                    include_space,
                    require_each_included_type,
                )
            }
        )

    def list_secret_version_ids(self, secret_id: str) -> str:
        secret = self.secrets[secret_id]

        version_list = []
        for version_id, version in secret.versions.items():
            version_list.append(
                {
                    "CreatedDate": int(time.time()),
                    "LastAccessedDate": int(time.time()),
                    "VersionId": version_id,
                    "VersionStages": version["version_stages"],
                }
            )

        return json.dumps(
            {
                "ARN": secret.secret_id,
                "Name": secret.name,
                "NextToken": "",
                "Versions": version_list,
            }
        )

    def list_secrets(
        self,
        filters: list[dict[str, Any]],
        max_results: int = MAX_RESULTS_DEFAULT,
        next_token: Optional[str] = None,
        include_planned_deletion: bool = False,
    ) -> tuple[list[dict[str, Any]], Optional[str]]:
        secret_list: list[dict[str, Any]] = []
        for secret in self.secrets.values():
            if hasattr(secret, "deleted_date"):
                if secret.deleted_date and not include_planned_deletion:
                    continue
            if _matches(secret, filters):
                secret_list.append(secret.to_dict())

        return self._get_secret_values_page_and_next_token(
            secret_list, max_results, next_token
        )

    def delete_secret(
        self,
        secret_id: str,
        recovery_window_in_days: Optional[int],
        force_delete_without_recovery: bool,
    ) -> tuple[str, str, float]:
        if recovery_window_in_days is not None and (
            recovery_window_in_days < 7 or recovery_window_in_days > 30
        ):
            raise InvalidParameterException(
                "An error occurred (InvalidParameterException) when calling the DeleteSecret operation: The \
                RecoveryWindowInDays value must be between 7 and 30 days (inclusive)."
            )

        if recovery_window_in_days and force_delete_without_recovery:
            raise InvalidParameterException(
                "An error occurred (InvalidParameterException) when calling the DeleteSecret operation: You can't \
                use ForceDeleteWithoutRecovery in conjunction with RecoveryWindowInDays."
            )

        if not self._is_valid_identifier(secret_id):
            if not force_delete_without_recovery:
                raise SecretNotFoundException()
            else:
                arn = SecretsManagerSecretIdentifier(
                    self.account_id, self.region_name, secret_id=secret_id
                ).generate()
                name = secret_id
                deletion_date = utcnow()
                return arn, name, self._unix_time_secs(deletion_date)
        else:
            secret = self.secrets[secret_id]
            if isinstance(secret, ReplicaSecret):
                raise OperationNotPermittedOnReplica
            if len(secret.replicas) > 0:
                replica_regions = ",".join([rep.region for rep in secret.replicas])
                msg = f"You can't delete secret {secret_id} that still has replica regions [{replica_regions}]"
                raise InvalidParameterException(msg)

            if secret.is_deleted() and not force_delete_without_recovery:
                raise InvalidRequestException(
                    "An error occurred (InvalidRequestException) when calling the DeleteSecret operation: You tried to \
                    perform the operation on a secret that's currently marked deleted."
                )

            deletion_date = utcnow()

            if force_delete_without_recovery:
                self.secrets.pop(secret_id)
            else:
                deletion_date += datetime.timedelta(days=recovery_window_in_days or 30)
                secret.delete(self._unix_time_secs(deletion_date))

            if not secret:
                raise SecretNotFoundException()

            arn = secret.arn
            name = secret.name

            return arn, name, self._unix_time_secs(deletion_date)

    def restore_secret(self, secret_id: str) -> tuple[str, str]:
        if not self._is_valid_identifier(secret_id):
            raise SecretNotFoundException()

        secret = self.secrets[secret_id]
        if isinstance(secret, ReplicaSecret):
            raise OperationNotPermittedOnReplica
        secret.restore()

        return secret.arn, secret.name

    def tag_resource(self, secret_id: str, tags: list[dict[str, str]]) -> None:
        if secret_id not in self.secrets:
            raise SecretNotFoundException()

        secret = self.secrets[secret_id]
        if isinstance(secret, ReplicaSecret):
            raise OperationNotPermittedOnReplica

        old_tags = {tag["Key"]: tag for tag in secret.tags or []}

        for tag in tags:
            old_tags[tag["Key"]] = tag

        secret.tags = list(old_tags.values())

    def untag_resource(self, secret_id: str, tag_keys: list[str]) -> None:
        if secret_id not in self.secrets:
            raise SecretNotFoundException()

        secret = self.secrets[secret_id]
        if isinstance(secret, ReplicaSecret):
            raise OperationNotPermittedOnReplica

        if secret.tags is None:
            return

        secret.tags = [tag for tag in secret.tags if tag["Key"] not in tag_keys]

    def update_secret_version_stage(
        self,
        secret_id: str,
        version_stage: str,
        remove_from_version_id: str,
        move_to_version_id: str,
    ) -> tuple[str, str]:
        if secret_id not in self.secrets:
            raise SecretNotFoundException()

        secret = self.secrets[secret_id]

        if remove_from_version_id:
            if remove_from_version_id not in secret.versions:
                raise InvalidParameterException(
                    f"Not a valid version: {remove_from_version_id}"
                )

            stages = secret.versions[remove_from_version_id]["version_stages"]
            if version_stage not in stages:
                raise InvalidParameterException(
                    f"Version stage {version_stage} not found in version {remove_from_version_id}"
                )

            stages.remove(version_stage)
        elif version_stage == "AWSCURRENT":
            current_version = [
                v
                for v in secret.versions
                if "AWSCURRENT" in secret.versions[v]["version_stages"]
            ][0]
            err = f"The parameter RemoveFromVersionId can't be empty. Staging label AWSCURRENT is currently attached to version {current_version}, so you must explicitly reference that version in RemoveFromVersionId."
            raise InvalidParameterException(err)

        if move_to_version_id:
            if move_to_version_id not in secret.versions:
                raise InvalidParameterException(
                    f"Not a valid version: {move_to_version_id}"
                )

            stages = secret.versions[move_to_version_id]["version_stages"]
            stages.append(version_stage)

        if version_stage == "AWSCURRENT":
            if remove_from_version_id:
                # Whenever you move AWSCURRENT, Secrets Manager automatically
                # moves the label AWSPREVIOUS to the version that AWSCURRENT
                # was removed from.
                for version in secret.versions:
                    if "AWSPREVIOUS" in secret.versions[version]["version_stages"]:
                        secret.versions[version]["version_stages"].remove("AWSPREVIOUS")
                secret.versions[remove_from_version_id]["version_stages"].append(
                    "AWSPREVIOUS"
                )

            if move_to_version_id:
                stages = secret.versions[move_to_version_id]["version_stages"]
                if "AWSPREVIOUS" in stages:
                    stages.remove("AWSPREVIOUS")

        return secret.arn, secret.name

    def put_resource_policy(self, secret_id: str, policy: str) -> tuple[str, str]:
        """
        The BlockPublicPolicy-parameter is not yet implemented
        """
        if not self._is_valid_identifier(secret_id):
            raise SecretNotFoundException()

        secret = self.secrets[secret_id]
        if isinstance(secret, ReplicaSecret):
            raise OperationNotPermittedOnReplica
        secret.policy = policy
        return secret.arn, secret.name

    def get_resource_policy(self, secret_id: str) -> str:
        if not self._is_valid_identifier(secret_id):
            raise SecretNotFoundException()

        secret = self.secrets[secret_id]
        if isinstance(secret, ReplicaSecret):
            raise OperationNotPermittedOnReplica
        resp = {
            "ARN": secret.arn,
            "Name": secret.name,
        }
        if secret.policy is not None:
            resp["ResourcePolicy"] = secret.policy
        return json.dumps(resp)

    def delete_resource_policy(self, secret_id: str) -> tuple[str, str]:
        if not self._is_valid_identifier(secret_id):
            raise SecretNotFoundException()

        secret = self.secrets[secret_id]
        if isinstance(secret, ReplicaSecret):
            raise OperationNotPermittedOnReplica
        secret.policy = None
        return secret.arn, secret.name

    def replicate_secret_to_regions(
        self,
        secret_id: str,
        replica_regions: list[dict[str, str]],
        force_overwrite: bool,
    ) -> tuple[str, list[dict[str, Any]]]:
        secret = self.describe_secret(secret_id)
        if isinstance(secret, ReplicaSecret):
            raise OperationNotPermittedOnReplica
        secret.replicas.extend(
            secret.create_replicas(replica_regions, force_overwrite=force_overwrite)
        )
        statuses = [replica.config for replica in secret.replicas]
        return secret_id, statuses

    def remove_regions_from_replication(
        self, secret_id: str, replica_regions: list[str]
    ) -> tuple[str, list[dict[str, str]]]:
        secret = self.describe_secret(secret_id)
        if isinstance(secret, ReplicaSecret):
            raise OperationNotPermittedOnReplica
        for replica in secret.replicas.copy():
            if replica.region in replica_regions:
                backend = secretsmanager_backends[self.account_id][replica.region]
                if replica.has_replica:
                    dict.pop(backend.secrets, replica.arn)
                secret.replicas.remove(replica)

        statuses = [replica.config for replica in secret.replicas]
        return secret_id, statuses

    def _get_secret_values_page_and_next_token(
        self,
        secret_list: list[dict[str, Any]],
        max_results: Optional[int],
        next_token: Optional[str],
    ) -> tuple[list[dict[str, Any]], Optional[str]]:
        starting_point = int(next_token or 0)
        ending_point = starting_point + int(max_results or MAX_RESULTS_DEFAULT)
        secret_page = secret_list[starting_point:ending_point]
        new_next_token = str(ending_point) if ending_point < len(secret_list) else None

        return secret_page, new_next_token


secretsmanager_backends = BackendDict(SecretsManagerBackend, "secretsmanager")
