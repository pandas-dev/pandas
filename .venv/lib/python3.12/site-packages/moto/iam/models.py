from __future__ import annotations

import json
import os
import re
import string
from collections.abc import Iterable
from datetime import date, datetime
from typing import Any, Optional, Union
from urllib import parse

from cryptography import x509
from cryptography.hazmat.backends import default_backend

from moto.core import DEFAULT_ACCOUNT_ID
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.exceptions import RESTError
from moto.core.utils import unix_time, utcnow
from moto.iam.policy_validation import (
    IAMPolicyDocumentValidator,
    IAMTrustPolicyDocumentValidator,
)
from moto.moto_api._internal import mock_random as random
from moto.settings import load_iam_aws_managed_policies
from moto.utilities.utils import (
    ARN_PARTITION_REGEX,
    PARTITION_NAMES,
    get_partition,
    md5_hash,
)

from ..utilities.tagging_service import TaggingService
from .aws_managed_policies import aws_managed_policies_data
from .exceptions import (
    DeleteConflictException,
    DuplicateTags,
    EntityAlreadyExists,
    InvalidInput,
    InvalidTagCharacters,
    LimitExceededException,
    MalformedCertificate,
    NoSuchEntity,
    NotFoundException,
    ReportNotPresentException,
    TagKeyTooBig,
    TagValueTooBig,
    TooManyTags,
    ValidationError,
)
from .utils import (
    generate_access_key_id_from_account_id,
    random_access_key,
    random_alphanumeric,
    random_policy_id,
    random_resource_id,
    random_role_id,
)

# Map to convert service names used in ServiceLinkedRoles
# The PascalCase should be used as part of the RoleName
SERVICE_NAME_CONVERSION = {
    "autoscaling": "AutoScaling",
    "application-autoscaling": "ApplicationAutoScaling",
    "elasticbeanstalk": "ElasticBeanstalk",
}


def get_account_id_from(access_key: str) -> str:
    # wrapped in a list() to avoid thread pooling problems (issue #5881)
    for account_id, account in list(iam_backends.items()):
        for partition in PARTITION_NAMES:
            if access_key in account[partition].access_keys:
                return account_id
    return DEFAULT_ACCOUNT_ID


def mark_account_as_visited(
    account_id: str, access_key: str, service: str, region: str
) -> None:
    partition = get_partition(region)
    account = iam_backends[account_id]
    if access_key in account[partition].access_keys:
        key = account[partition].access_keys[access_key]
        key.last_used = AccessKeyLastUsed(
            timestamp=utcnow(), service=service, region=region
        )
        if key.role_arn:
            try:
                role = account[partition].get_role_by_arn(key.role_arn)
                role.last_used = utcnow()
            except NotFoundException:
                # User assumes a non-existing role
                pass
    else:
        # User provided access credentials unknown to us
        pass


def _serialize_version_datetime(value: Any) -> str:
    if isinstance(value, date):
        return value.strftime("%Y-%m-%d")
    raise TypeError("Unable to serialize value.")


LIMIT_KEYS_PER_USER = 2


class MFADevice:
    """MFA Device class."""

    def __init__(
        self, serial_number: str, authentication_code_1: str, authentication_code_2: str
    ):
        self.enable_date = utcnow()
        self.serial_number = serial_number
        self.authentication_code_1 = authentication_code_1
        self.authentication_code_2 = authentication_code_2


class VirtualMfaDevice:
    def __init__(self, account_id: str, region_name: str, device_name: str):
        self.serial_number = (
            f"arn:{get_partition(region_name)}:iam::{account_id}:mfa{device_name}"
        )

        random_base32_string = "".join(
            random.choice(string.ascii_uppercase + "234567") for _ in range(64)
        )
        self.base32_string_seed = random_base32_string
        self.qr_code_png = os.urandom(64)  # this would be a generated PNG

        self.enable_date: Optional[datetime] = None
        self.user_attribute: Optional[dict[str, Any]] = None
        self.user: Optional[User] = None


class Policy(CloudFormationModel):
    # Note: This class does not implement the CloudFormation support for AWS::IAM::Policy, as that CF resource
    #  is for creating *inline* policies.  That is done in class InlinePolicy.

    is_attachable = False

    def __init__(
        self,
        name: str,
        account_id: str,
        region: str,
        default_version_id: Optional[str] = None,
        description: Optional[str] = None,
        document: Optional[str] = None,
        path: Optional[str] = None,
        create_date: Optional[datetime] = None,
        update_date: Optional[datetime] = None,
        tags: Optional[dict[str, dict[str, str]]] = None,
    ):
        self.name = name
        self.account_id = account_id
        self.attachment_count = 0
        self.description = description or ""
        self.id = random_policy_id()
        self.path = path or "/"
        self.tags = tags or {}
        self.partition = get_partition(region)

        if default_version_id:
            self.default_version_id = default_version_id
            self.next_version_num = int(default_version_id.lstrip("v")) + 1
        else:
            self.default_version_id = "v1"
            self.next_version_num = 2
        self.versions = [
            PolicyVersion(
                self.arn,  # type: ignore[attr-defined]
                document,
                True,
                self.default_version_id,
                update_date,
            )
        ]

        self.create_date = create_date or utcnow()
        self.update_date = update_date or utcnow()

    def update_default_version(self, new_default_version_id: str) -> None:
        for version in self.versions:
            if version.version_id == new_default_version_id:
                version.is_default = True
            if version.version_id == self.default_version_id:
                version.is_default = False
        self.default_version_id = new_default_version_id

    @property
    def policy_version_list(self) -> list[dict[str, bool | str | datetime | None]]:
        return [
            {
                "Document": parse.quote(version.document),
                "VersionId": version.version_id,
                "IsDefaultVersion": version.is_default_version,
                "CreateDate": version.create_date,
            }
            for version in self.versions
        ]

    def get_tags(self) -> list[dict[str, str]]:
        return [self.tags[tag] for tag in self.tags]

    @property
    def tag_list_type(self) -> list[dict[str, str]]:
        return self.get_tags()


class SAMLProvider(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        saml_metadata_document: Optional[str] = None,
    ):
        self.account_id = account_id
        self.create_date = utcnow()
        self.region_name = region_name
        self.name = name
        self.saml_metadata_document = saml_metadata_document
        self.valid_until = datetime(2032, 5, 9, 16, 27, 11)

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.region_name)}:iam::{self.account_id}:saml-provider/{self.name}"


class OpenIDConnectProvider(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        url: str,
        thumbprint_list: list[str],
        client_id_list: list[str],
        tags: dict[str, dict[str, str]],
    ):
        self._errors: list[str] = []
        self._validate(url, thumbprint_list, client_id_list)

        self.account_id = account_id
        self.region_name = region_name
        parsed_url = parse.urlparse(url)
        self.url = parsed_url.netloc + parsed_url.path
        self.thumbprint_list = thumbprint_list
        self.client_id_list = client_id_list
        self.create_date = utcnow()
        self.tags = tags or {}

    @property
    def tag_list_type(self) -> list[dict[str, str]]:
        tags = []
        for key in self.tags:
            tags.append(self.tags[key])
        return tags

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.region_name)}:iam::{self.account_id}:oidc-provider/{self.url}"

    def _validate(
        self, url: str, thumbprint_list: list[str], client_id_list: list[str]
    ) -> None:
        if any(len(client_id) > 255 for client_id in client_id_list):
            self._errors.append(
                self._format_error(
                    key="clientIDList",
                    value=client_id_list,
                    constraint="Member must satisfy constraint: "
                    "[Member must have length less than or equal to 255, "
                    "Member must have length greater than or equal to 1]",
                )
            )

        if any(len(thumbprint) > 40 for thumbprint in thumbprint_list):
            self._errors.append(
                self._format_error(
                    key="thumbprintList",
                    value=thumbprint_list,
                    constraint="Member must satisfy constraint: "
                    "[Member must have length less than or equal to 40, "
                    "Member must have length greater than or equal to 40]",
                )
            )

        if len(url) > 255:
            self._errors.append(
                self._format_error(
                    key="url",
                    value=url,
                    constraint="Member must have length less than or equal to 255",
                )
            )

        self._raise_errors()

        parsed_url = parse.urlparse(url)
        if not parsed_url.scheme or not parsed_url.netloc:
            raise ValidationError("Invalid Open ID Connect Provider URL")

        if len(thumbprint_list) > 5:
            raise InvalidInput("Thumbprint list must contain fewer than 5 entries.")

        if len(client_id_list) > 100:
            raise LimitExceededException(
                "Cannot exceed quota for ClientIdsPerOpenIdConnectProvider: 100"
            )

    def _format_error(self, key: str, value: Any, constraint: str) -> str:
        return f'Value "{value}" at "{key}" failed to satisfy constraint: {constraint}'

    def _raise_errors(self) -> None:
        if self._errors:
            count = len(self._errors)
            plural = "s" if len(self._errors) > 1 else ""
            errors = "; ".join(self._errors)
            self._errors = []  # reset collected errors

            raise ValidationError(
                f"{count} validation error{plural} detected: {errors}"
            )

    def get_tags(self) -> list[dict[str, str]]:
        return [self.tags[tag] for tag in self.tags]


class PolicyVersion:
    def __init__(
        self,
        policy_arn: str,
        document: Optional[str],
        is_default: bool = False,
        version_id: str = "v1",
        create_date: Optional[datetime] = None,
    ):
        self.policy_arn = policy_arn
        self.document = document or ""
        self.is_default = is_default
        self.version_id = version_id

        self.create_date = create_date or utcnow()

    @property
    def is_default_version(self) -> bool:
        return self.is_default


class ManagedPolicy(Policy, CloudFormationModel):
    """Managed policy."""

    @property
    def backend(self) -> IAMBackend:
        return iam_backends[self.account_id][self.partition]

    is_attachable = True

    def attach_to(self, obj: Union[Role, Group, User]) -> None:
        self.attachment_count += 1
        obj.managed_policies[self.arn] = self

    def detach_from(self, obj: Union[Role, Group, User]) -> None:
        self.attachment_count -= 1
        del obj.managed_policies[self.arn]

    @property
    def arn(self) -> str:
        return (
            f"arn:{self.partition}:iam::{self.account_id}:policy{self.path}{self.name}"
        )

    def to_config_dict(self) -> dict[str, Any]:
        return {
            "version": "1.3",
            "configurationItemCaptureTime": str(self.create_date),
            "configurationItemStatus": "OK",
            "configurationStateId": str(int(unix_time())),
            "arn": f"arn:{self.partition}:iam::{self.account_id}:policy/{self.name}",
            "resourceType": "AWS::IAM::Policy",
            "resourceId": self.id,
            "resourceName": self.name,
            "awsRegion": "global",
            "availabilityZone": "Not Applicable",
            "resourceCreationTime": str(self.create_date),
            "tags": self.tags,
            "configuration": {
                "policyName": self.name,
                "policyId": self.id,
                "arn": f"arn:{self.partition}:iam::{self.account_id}:policy/{self.name}",
                "path": self.path,
                "defaultVersionId": self.default_version_id,
                "attachmentCount": self.attachment_count,
                "permissionsBoundaryUsageCount": 0,
                "isAttachable": ManagedPolicy.is_attachable,
                "description": self.description,
                "createDate": str(self.create_date.isoformat()),
                "updateDate": str(self.create_date.isoformat()),
                "tags": [
                    {"key": key, "value": self.tags[key]["Value"]} for key in self.tags
                ],
                "policyVersionList": [
                    {
                        "document": parse.quote(version.document),
                        "versionId": version.version_id,
                        "isDefaultVersion": version.is_default,
                        "createDate": str(version.create_date),
                    }
                    for version in self.versions
                ],
            },
            "supplementaryConfiguration": {},
        }

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""  # Resource never gets named after by template PolicyName!

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::IAM::ManagedPolicy"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> ManagedPolicy:
        properties = cloudformation_json.get("Properties", {})
        policy_document = json.dumps(properties.get("PolicyDocument"))
        name = properties.get("ManagedPolicyName", resource_name)
        description = properties.get("Description")
        path = properties.get("Path")
        group_names = properties.get("Groups", [])
        user_names = properties.get("Users", [])
        role_names = properties.get("Roles", [])
        tags = properties.get("Tags", {})

        partition = get_partition(region_name)
        policy = iam_backends[account_id][partition].create_policy(
            description=description,
            path=path,
            policy_document=policy_document,
            policy_name=name,
            tags=tags,
        )
        for group_name in group_names:
            iam_backends[account_id][partition].attach_group_policy(
                group_name=group_name, policy_arn=policy.arn
            )
        for user_name in user_names:
            iam_backends[account_id][partition].attach_user_policy(
                user_name=user_name, policy_arn=policy.arn
            )
        for role_name in role_names:
            iam_backends[account_id][partition].attach_role_policy(
                role_name=role_name, policy_arn=policy.arn
            )
        return policy

    def __eq__(self, other: Any) -> bool:
        return self.arn == other.arn

    def __hash__(self) -> int:
        return self.arn.__hash__()

    @property
    def physical_resource_id(self) -> str:
        return self.arn


class AWSManagedPolicy(ManagedPolicy):
    """AWS-managed policy."""

    @classmethod
    def from_data(  # type: ignore[misc]
        cls, name: str, account_id: str, region_name: str, data: dict[str, Any]
    ) -> AWSManagedPolicy:
        return cls(
            name,
            account_id=account_id,
            region=region_name,
            default_version_id=data.get("DefaultVersionId"),
            path=data.get("Path"),
            document=json.dumps(data.get("Document")),
            create_date=datetime.fromisoformat(data["CreateDate"]),
            update_date=datetime.fromisoformat(data["UpdateDate"]),
        )

    @property
    def arn(self) -> str:
        return f"arn:{self.partition}:iam::aws:policy{self.path}{self.name}"


class InlinePolicy(CloudFormationModel):
    # Represents an Inline Policy created by CloudFormation
    def __init__(
        self,
        resource_name: str,
        policy_name: str,
        policy_document: str,
        group_names: list[str],
        role_names: list[str],
        user_names: list[str],
    ):
        self.name = resource_name
        self.policy_name = policy_name
        self.policy_document = policy_document
        self.group_names = group_names
        self.role_names = role_names
        self.user_names = user_names
        self.update(policy_name, policy_document, group_names, role_names, user_names)

    def update(
        self,
        policy_name: str,
        policy_document: str,
        group_names: list[str],
        role_names: list[str],
        user_names: list[str],
    ) -> None:
        self.policy_name = policy_name
        self.policy_document = (
            json.dumps(policy_document)
            if isinstance(policy_document, dict)
            else policy_document
        )
        self.group_names = group_names
        self.role_names = role_names
        self.user_names = user_names

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""  # Resource never gets named after by template PolicyName!

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::IAM::Policy"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> InlinePolicy:
        properties = cloudformation_json.get("Properties", {})
        policy_document = properties.get("PolicyDocument")
        policy_name = properties.get("PolicyName")
        user_names = properties.get("Users")
        role_names = properties.get("Roles")
        group_names = properties.get("Groups")

        return iam_backends[account_id][
            get_partition(region_name)
        ].create_inline_policy(
            resource_name,
            policy_name,
            policy_document,
            group_names,
            role_names,
            user_names,
        )

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: InlinePolicy,
        new_resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> InlinePolicy:
        properties = cloudformation_json["Properties"]

        if cls.is_replacement_update(properties):
            resource_name_property = cls.cloudformation_name_type()
            if resource_name_property not in properties:
                properties[resource_name_property] = new_resource_name
            new_resource = cls.create_from_cloudformation_json(
                properties[resource_name_property],
                cloudformation_json,
                account_id,
                region_name,
            )
            properties[resource_name_property] = original_resource.name
            cls.delete_from_cloudformation_json(
                original_resource.name, cloudformation_json, account_id, region_name
            )
            return new_resource

        else:  # No Interruption
            properties = cloudformation_json.get("Properties", {})
            policy_document = properties.get("PolicyDocument")
            policy_name = properties.get("PolicyName", original_resource.name)
            user_names = properties.get("Users")
            role_names = properties.get("Roles")
            group_names = properties.get("Groups")

            return iam_backends[account_id][
                get_partition(region_name)
            ].update_inline_policy(
                original_resource.name,
                policy_name,
                policy_document,
                group_names,
                role_names,
                user_names,
            )

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> None:
        iam_backends[account_id][get_partition(region_name)].delete_inline_policy(
            resource_name
        )

    @staticmethod
    def is_replacement_update(properties: list[str]) -> bool:
        properties_requiring_replacement_update: list[str] = []
        return any(
            property_requiring_replacement in properties
            for property_requiring_replacement in properties_requiring_replacement_update
        )

    @property
    def physical_resource_id(self) -> str:
        return self.name

    def apply_policy(self, backend: IAMBackend) -> None:
        if self.user_names:
            for user_name in self.user_names:
                backend.put_user_policy(
                    user_name, self.policy_name, self.policy_document
                )
        if self.role_names:
            for role_name in self.role_names:
                backend.put_role_policy(
                    role_name, self.policy_name, self.policy_document
                )
        if self.group_names:
            for group_name in self.group_names:
                backend.put_group_policy(
                    group_name, self.policy_name, self.policy_document
                )

    def unapply_policy(self, backend: IAMBackend) -> None:
        if self.user_names:
            for user_name in self.user_names:
                backend.delete_user_policy(user_name, self.policy_name)
        if self.role_names:
            for role_name in self.role_names:
                backend.delete_role_policy(role_name, self.policy_name)
        if self.group_names:
            for group_name in self.group_names:
                backend.delete_group_policy(group_name, self.policy_name)


class Role(CloudFormationModel):
    def __init__(
        self,
        account_id: str,
        partition: str,
        role_id: str,
        name: str,
        assume_role_policy_document: str,
        path: str,
        permissions_boundary: Optional[str],
        description: str,
        tags: list[dict[str, str]],
        max_session_duration: Optional[str],
        linked_service: Optional[str] = None,
    ):
        self.account_id = account_id
        self.partition = partition
        self.id = role_id
        self.name = name
        self.assume_role_policy_document = assume_role_policy_document
        self.path = path or "/"
        self.policies: dict[str, str] = {}
        self.managed_policies: dict[str, ManagedPolicy] = {}
        self.create_date = utcnow()
        self.tags = tags
        # last_used should be treated as part of the public API
        # https://github.com/getmoto/moto/issues/6859
        self.last_used: Optional[datetime] = None
        self.last_used_region = None
        self.description = description
        self.permissions_boundary_arn: Optional[str] = permissions_boundary
        self.max_session_duration = max_session_duration
        self._linked_service = linked_service

    @property
    def attached_managed_policies(self) -> list[dict[str, str]]:
        managed_policies = []
        for policy in self.managed_policies.values():
            managed_policies.append(
                {
                    "PolicyName": policy.name,
                    "PolicyArn": policy.arn,
                }
            )
        return managed_policies

    @property
    def instance_profile_list(self) -> list[InstanceProfile]:
        backend = iam_backends[self.account_id][self.partition]
        instance_profiles = []
        for instance_profile in backend.instance_profiles.values():
            if self in instance_profile.roles:
                instance_profiles.append(instance_profile)
        return instance_profiles

    @property
    def role_policy_list(self) -> list[dict[str, str]]:
        policy_list = [
            {"PolicyName": name, "PolicyDocument": parse.quote(policy)}
            for name, policy in self.policies.items()
        ]
        return policy_list

    @property
    def permissions_boundary(self) -> Optional[dict[str, str]]:
        if self.permissions_boundary_arn:
            return {
                "PermissionsBoundaryType": "PermissionsBoundaryPolicy",
                "PermissionsBoundaryArn": self.permissions_boundary_arn,
            }
        return None

    @property
    def role_last_used(self) -> dict[str, Optional[str | datetime]]:
        """Returns a dict with the last used information of the role."""
        return {
            "LastUsedDate": self.last_used,
            "Region": self.last_used_region,
        }

    @staticmethod
    def cloudformation_name_type() -> str:
        return "RoleName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-iam-role.html
        return "AWS::IAM::Role"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> Role:
        properties = cloudformation_json["Properties"]
        role_name = properties.get("RoleName", resource_name)

        assume_role_policy_document = properties["AssumeRolePolicyDocument"]
        if not isinstance(assume_role_policy_document, str):
            assume_role_policy_document = json.dumps(
                assume_role_policy_document, default=_serialize_version_datetime
            )

        iam_backend = iam_backends[account_id][get_partition(region_name)]
        role = iam_backend.create_role(
            role_name=role_name,
            assume_role_policy_document=assume_role_policy_document,
            path=properties.get("Path", "/"),
            permissions_boundary=properties.get("PermissionsBoundary", ""),
            description=properties.get("Description", ""),
            tags=properties.get("Tags", {}),
            max_session_duration=properties.get("MaxSessionDuration", 3600),
        )

        policies = properties.get("Policies", [])
        for policy in policies:
            policy_name = policy["PolicyName"]
            policy_json = policy["PolicyDocument"]
            role.put_policy(policy_name, policy_json)

        return role

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Role,
        new_resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> Role:
        iam_backend = iam_backends[account_id][get_partition(region_name)]
        properties = cloudformation_json["Properties"]
        role_name = properties.get("RoleName", original_resource.name)
        role_description = properties.get("Description", "")
        max_session_duration = properties.get("RoleDescription", 3600)
        return iam_backend.update_role(
            role_name, role_description, max_session_duration
        )

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> None:
        backend: IAMBackend = iam_backends[account_id][get_partition(region_name)]
        for profile in backend.instance_profiles.values():
            profile.delete_role(role_name=resource_name)

        for role in backend.roles.values():
            if role.name == resource_name:
                for arn in list(role.policies.keys()):
                    role.delete_policy(arn)

        # We may have inline policies attached to this role - make sure that they don't reference a deleted role
        # If the reference still exists, it's impossible to delete the policy, because on deletion
        # it will try to detach itself from a role that no longer exists
        for policy in backend.inline_policies.values():
            policy.role_names = [n for n in policy.role_names if n != resource_name]

        backend.delete_role(resource_name)

    @property
    def arn(self) -> str:
        if self._linked_service:
            return f"arn:{self.partition}:iam::{self.account_id}:role/aws-service-role/{self._linked_service}/{self.name}"
        return f"arn:{self.partition}:iam::{self.account_id}:role{self.path}{self.name}"

    def to_config_dict(self) -> dict[str, Any]:
        _managed_policies = []
        for key in self.managed_policies.keys():
            _managed_policies.append(
                {
                    "policyArn": key,
                    "policyName": iam_backends[self.account_id][self.partition]
                    .managed_policies[key]
                    .name,
                }
            )

        _role_policy_list = []
        for key, value in self.policies.items():
            _role_policy_list.append(
                {"policyName": key, "policyDocument": parse.quote(value)}
            )

        _instance_profiles = []
        backend = iam_backends[self.account_id][self.partition]
        for instance_profile in backend.instance_profiles.values():
            for _ in instance_profile.roles:
                _instance_profiles.append(instance_profile.to_embedded_config_dict())
                break

        config_dict = {
            "version": "1.3",
            "configurationItemCaptureTime": str(self.create_date),
            "configurationItemStatus": "ResourceDiscovered",
            "configurationStateId": str(int(unix_time())),
            "arn": f"arn:{self.partition}:iam::{self.account_id}:role/{self.name}",
            "resourceType": "AWS::IAM::Role",
            "resourceId": self.name,
            "resourceName": self.name,
            "awsRegion": "global",
            "availabilityZone": "Not Applicable",
            "resourceCreationTime": str(self.create_date),
            "relatedEvents": [],
            "relationships": [],
            "tags": self.tags,
            "configuration": {
                "path": self.path,
                "roleName": self.name,
                "roleId": self.id,
                "arn": f"arn:{self.partition}:iam::{self.account_id}:role/{self.name}",
                "assumeRolePolicyDocument": parse.quote(
                    self.assume_role_policy_document
                )
                if self.assume_role_policy_document
                else None,
                "instanceProfileList": _instance_profiles,
                "rolePolicyList": _role_policy_list,
                "createDate": self.create_date.isoformat(),
                "attachedManagedPolicies": _managed_policies,
                "permissionsBoundary": self.permissions_boundary_arn,
                "tags": [
                    {"key": tag["Key"], "value": tag["Value"]} for tag in self.tags
                ],
                "roleLastUsed": None,
            },
            "supplementaryConfiguration": {},
        }
        return config_dict

    def put_policy(self, policy_name: str, policy_json: str) -> None:
        self.policies[policy_name] = policy_json

    def delete_policy(self, policy_name: str) -> None:
        try:
            del self.policies[policy_name]
        except KeyError:
            raise NotFoundException(
                f"The role policy with name {policy_name} cannot be found."
            )

    @property
    def physical_resource_id(self) -> str:
        return self.name

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn", "RoleId"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return self.arn
        if attribute_name == "RoleId":
            return self.id
        raise UnformattedGetAttTemplateException()

    @property
    def description_escaped(self) -> str:
        import html

        return html.escape(self.description or "")


class InstanceProfile(CloudFormationModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        instance_profile_id: str,
        name: str,
        path: str,
        roles: list[Role],
        tags: Optional[list[dict[str, str]]] = None,
    ):
        self.id = instance_profile_id
        self.account_id = account_id
        self.partition = get_partition(region_name)
        self.name = name
        self.path = path or "/"
        self.roles = roles if roles else []
        self.create_date = utcnow()
        self.tags = tags or []

    @staticmethod
    def cloudformation_name_type() -> str:
        return "InstanceProfileName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-iam-instanceprofile.html
        return "AWS::IAM::InstanceProfile"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> InstanceProfile:
        properties = cloudformation_json["Properties"]

        role_names = properties["Roles"]
        return iam_backends[account_id][
            get_partition(region_name)
        ].create_instance_profile(
            name=resource_name,
            path=properties.get("Path", "/"),
            role_names=role_names,
        )

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> None:
        iam_backends[account_id][get_partition(region_name)].delete_instance_profile(
            resource_name, ignore_attached_roles=True
        )

    def delete_role(self, role_name: str) -> None:
        self.roles = [role for role in self.roles if role.name != role_name]

    @property
    def arn(self) -> str:
        return f"arn:{self.partition}:iam::{self.account_id}:instance-profile{self.path}{self.name}"

    @property
    def physical_resource_id(self) -> str:
        return self.name

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return self.arn
        raise UnformattedGetAttTemplateException()

    def to_embedded_config_dict(self) -> dict[str, Any]:
        # Instance Profiles aren't a config item itself, but they are returned in IAM roles with
        # a "config like" json structure. It's also different than Role.to_config_dict()
        roles = []
        for role in self.roles:
            roles.append(
                {
                    "path": role.path,
                    "roleName": role.name,
                    "roleId": role.id,
                    "arn": f"arn:{self.partition}:iam::{self.account_id}:role/{role.name}",
                    "createDate": str(role.create_date),
                    "assumeRolePolicyDocument": parse.quote(
                        role.assume_role_policy_document
                    ),
                    "description": role.description,
                    "maxSessionDuration": None,
                    "permissionsBoundary": role.permissions_boundary_arn,
                    "tags": role.tags,
                    "roleLastUsed": None,
                }
            )

        return {
            "path": self.path,
            "instanceProfileName": self.name,
            "instanceProfileId": self.id,
            "arn": f"arn:{self.partition}:iam::{self.account_id}:instance-profile/{role.name}",
            "createDate": str(self.create_date),
            "roles": roles,
        }


class Certificate(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        cert_name: str,
        cert_body: str,
        private_key: str,
        cert_chain: Optional[str] = None,
        path: Optional[str] = None,
    ):
        self.account_id = account_id
        self.partition = get_partition(region_name)
        self.cert_name = cert_name
        if cert_body:
            cert_body = cert_body.rstrip()
        self.cert_body = cert_body
        self.private_key = private_key
        self.path = path if path else "/"
        self.cert_chain = cert_chain

    @property
    def physical_resource_id(self) -> str:
        return self.cert_name

    @property
    def arn(self) -> str:
        return f"arn:{self.partition}:iam::{self.account_id}:server-certificate{self.path}{self.cert_name}"


class SigningCertificate(BaseModel):
    def __init__(self, certificate_id: str, user_name: str, body: str):
        self.id = certificate_id
        self.user_name = user_name
        self.body = body
        self.upload_date = utcnow()
        self.status = "Active"

    @property
    def certificate_id(self) -> str:
        return self.id

    @property
    def certificate_body(self) -> str:
        return self.body


class AccessKeyLastUsed:
    def __init__(self, timestamp: datetime, service: str, region: str):
        self._timestamp = timestamp
        self.service = service
        self.region = region

    @property
    def timestamp(self) -> datetime:
        return self._timestamp

    def strftime(self, date_format: str) -> str:
        return self._timestamp.strftime(date_format)


class AccessKey(CloudFormationModel):
    def __init__(
        self,
        user_name: Optional[str],
        prefix: str,
        account_id: str,
        status: str = "Active",
    ):
        self.user_name = user_name
        self.access_key_id = generate_access_key_id_from_account_id(
            account_id, prefix=prefix, total_length=20
        )
        self.secret_access_key = random_alphanumeric(40)
        self.status = status
        self.create_date = utcnow()

        # Some users will set this field manually
        # And they will be setting this value to a `datetime`
        # https://github.com/getmoto/moto/issues/5927#issuecomment-1738188283
        #
        # The `to_csv` method calls `last_used.strptime`, which currently works on both AccessKeyLastUsed and datetime
        # In the next major release we should communicate that this only accepts AccessKeyLastUsed
        # (And rework to_csv accordingly)
        self.last_used: Optional[AccessKeyLastUsed] = None
        self.role_arn: Optional[str] = None

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["SecretAccessKey"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "SecretAccessKey":
            return self.secret_access_key
        raise UnformattedGetAttTemplateException()

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""  # Resource never gets named after by template PolicyName!

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::IAM::AccessKey"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> AccessKey:
        properties = cloudformation_json.get("Properties", {})
        user_name = properties.get("UserName")
        status = properties.get("Status", "Active")

        return iam_backends[account_id][get_partition(region_name)].create_access_key(
            user_name, status=status
        )

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: AccessKey,
        new_resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> AccessKey:
        properties = cloudformation_json["Properties"]

        if cls.is_replacement_update(properties):
            new_resource = cls.create_from_cloudformation_json(
                new_resource_name, cloudformation_json, account_id, region_name
            )
            cls.delete_from_cloudformation_json(
                original_resource.physical_resource_id,
                cloudformation_json,
                account_id,
                region_name,
            )
            return new_resource

        else:  # No Interruption
            properties = cloudformation_json.get("Properties", {})
            status = properties.get("Status")
            return iam_backends[account_id][
                get_partition(region_name)
            ].update_access_key(
                original_resource.user_name,  # type: ignore[arg-type]
                original_resource.access_key_id,
                status,
            )

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> None:
        iam_backends[account_id][get_partition(region_name)].delete_access_key_by_name(
            resource_name
        )

    @staticmethod
    def is_replacement_update(properties: list[str]) -> bool:
        properties_requiring_replacement_update = ["Serial", "UserName"]
        return any(
            property_requiring_replacement in properties
            for property_requiring_replacement in properties_requiring_replacement_update
        )

    @property
    def physical_resource_id(self) -> str:
        return self.access_key_id


class SshPublicKey(BaseModel):
    def __init__(self, user_name: str, ssh_public_key_body: str):
        self.user_name = user_name
        self.ssh_public_key_body = ssh_public_key_body
        self.ssh_public_key_id = "APKA" + random_access_key()
        self.fingerprint = md5_hash(ssh_public_key_body.encode()).hexdigest()
        self.status = "Active"
        self.upload_date = utcnow()


class Group(BaseModel):
    def __init__(self, account_id: str, region_name: str, name: str, path: str = "/"):
        self.account_id = account_id
        self.partition = get_partition(region_name)
        self.name = name
        self.id = random_resource_id()
        self.path = path
        self.create_date = utcnow()

        self.users: list[User] = []
        self.managed_policies: dict[str, ManagedPolicy] = {}
        self.policies: dict[str, str] = {}

    @property
    def group_policy_list(self) -> list[dict[str, str]]:
        policy_list = [
            {"PolicyName": name, "PolicyDocument": parse.quote(policy)}
            for name, policy in self.policies.items()
        ]
        return policy_list

    @property
    def attached_managed_policies(self) -> list[dict[str, str]]:
        attached_policies = []
        for policy in self.managed_policies.values():
            attached_policies.append(
                {
                    "PolicyName": policy.name,
                    "PolicyArn": policy.arn,
                }
            )
        return attached_policies

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn"]

    def get_cfn_attribute(self, attribute_name: str) -> None:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            raise NotImplementedError('"Fn::GetAtt" : [ "{0}" , "Arn" ]"')
        raise UnformattedGetAttTemplateException()

    @property
    def arn(self) -> str:
        if self.path == "/":
            return f"arn:{self.partition}:iam::{self.account_id}:group/{self.name}"
        else:
            # The path must by definition end and start with a forward slash. So we don't have to add more slashes to the ARN
            return f"arn:{self.partition}:iam::{self.account_id}:group{self.path}{self.name}"

    def get_policy(self, policy_name: str) -> dict[str, str]:
        try:
            policy_json = self.policies[policy_name]
        except KeyError:
            raise NotFoundException(f"Policy {policy_name} not found")

        return {
            "policy_name": policy_name,
            "policy_document": policy_json,
            "group_name": self.name,
        }

    def put_policy(self, policy_name: str, policy_json: str) -> None:
        self.policies[policy_name] = policy_json

    def list_policies(self) -> list[str]:
        return list(self.policies.keys())

    def delete_policy(self, policy_name: str) -> None:
        if policy_name not in self.policies:
            raise NotFoundException(f"Policy {policy_name} not found")

        del self.policies[policy_name]


class User(CloudFormationModel):
    def __init__(
        self, account_id: str, region_name: str, name: str, path: Optional[str] = None
    ):
        self.account_id = account_id
        self.region_name = region_name
        self.name = name
        self.id = random_resource_id()
        self.path = path if path else "/"
        self.create_date = utcnow()
        self.mfa_devices: dict[str, MFADevice] = {}
        self.policies: dict[str, str] = {}
        self.managed_policies: dict[str, ManagedPolicy] = {}
        self.access_keys: list[AccessKey] = []
        self.ssh_public_keys: list[SshPublicKey] = []
        self.password: Optional[str] = None
        # last_used should be treated as part of the public API
        # https://github.com/getmoto/moto/issues/5927
        self.password_last_used = None
        self.password_reset_required = False
        self.signing_certificates: dict[str, SigningCertificate] = {}

    @property
    def attached_managed_policies(self) -> list[dict[str, str]]:
        managed_policies = []
        for policy in self.managed_policies.values():
            managed_policies.append(
                {
                    "PolicyName": policy.name,
                    "PolicyArn": policy.arn,
                }
            )
        return managed_policies

    @property
    def user_policy_list(self) -> list[dict[str, str]]:
        policy_list = [
            {"PolicyName": name, "PolicyDocument": parse.quote(policy)}
            for name, policy in self.policies.items()
        ]
        return policy_list

    @property
    def group_list(self) -> list[str]:
        backend = iam_backends[self.account_id][get_partition(self.region_name)]
        groups = backend.get_groups_for_user(self.name)
        return [group.name for group in groups]

    @property
    def tags(self) -> Optional[list[dict[str, str]]]:
        backend = iam_backends[self.account_id][get_partition(self.region_name)]
        tags = backend.tagger.list_tags_for_resource(self.arn)["Tags"]
        return tags if tags else None

    @property
    def arn(self) -> str:
        partition = get_partition(self.region_name)
        return f"arn:{partition}:iam::{self.account_id}:user{self.path}{self.name}"

    def get_policy(self, policy_name: str) -> dict[str, str]:
        try:
            policy_json = self.policies[policy_name]
        except KeyError:
            raise NotFoundException(f"Policy {policy_name} not found")

        return {
            "policy_name": policy_name,
            "policy_document": policy_json,
            "user_name": self.name,
        }

    def put_policy(self, policy_name: str, policy_json: str) -> None:
        self.policies[policy_name] = policy_json

    def deactivate_mfa_device(self, serial_number: str) -> None:
        self.mfa_devices.pop(serial_number)

    def delete_policy(self, policy_name: str) -> None:
        if policy_name not in self.policies:
            raise NotFoundException(f"Policy {policy_name} not found")

        del self.policies[policy_name]

    def create_access_key(self, prefix: str, status: str = "Active") -> AccessKey:
        access_key = AccessKey(
            self.name, prefix=prefix, status=status, account_id=self.account_id
        )
        self.access_keys.append(access_key)
        return access_key

    def enable_mfa_device(
        self, serial_number: str, authentication_code_1: str, authentication_code_2: str
    ) -> None:
        self.mfa_devices[serial_number] = MFADevice(
            serial_number, authentication_code_1, authentication_code_2
        )

    def get_all_access_keys(self) -> list[AccessKey]:
        return self.access_keys

    def delete_access_key(self, access_key_id: str) -> None:
        key = self.get_access_key_by_id(access_key_id)
        self.access_keys.remove(key)

    def update_access_key(
        self, access_key_id: str, status: Optional[str] = None
    ) -> AccessKey:
        key = self.get_access_key_by_id(access_key_id)
        if status is not None:
            key.status = status
        return key

    def get_access_key_by_id(self, access_key_id: str) -> AccessKey:
        for key in self.access_keys:
            if key.access_key_id == access_key_id:
                return key

        raise NotFoundException(
            f"The Access Key with id {access_key_id} cannot be found"
        )

    def has_access_key(self, access_key_id: str) -> bool:
        return any(
            access_key
            for access_key in self.access_keys
            if access_key.access_key_id == access_key_id
        )

    def upload_ssh_public_key(self, ssh_public_key_body: str) -> SshPublicKey:
        pubkey = SshPublicKey(self.name, ssh_public_key_body)
        self.ssh_public_keys.append(pubkey)
        return pubkey

    def get_ssh_public_key(self, ssh_public_key_id: str) -> SshPublicKey:
        for key in self.ssh_public_keys:
            if key.ssh_public_key_id == ssh_public_key_id:
                return key

        raise NotFoundException(
            f"The SSH Public Key with id {ssh_public_key_id} cannot be found"
        )

    def get_all_ssh_public_keys(self) -> list[SshPublicKey]:
        return self.ssh_public_keys

    def update_ssh_public_key(self, ssh_public_key_id: str, status: str) -> None:
        key = self.get_ssh_public_key(ssh_public_key_id)
        key.status = status

    def delete_ssh_public_key(self, ssh_public_key_id: str) -> None:
        key = self.get_ssh_public_key(ssh_public_key_id)
        self.ssh_public_keys.remove(key)

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return self.arn
        raise UnformattedGetAttTemplateException()

    def to_csv(self) -> str:
        date_format = "%Y-%m-%dT%H:%M:%S+00:00"
        date_created = self.create_date
        # aagrawal,arn:aws:iam::509284790694:user/aagrawal,2014-09-01T22:28:48+00:00,true,2014-11-12T23:36:49+00:00,2014-09-03T18:59:00+00:00,N/A,false,true,2014-09-01T22:28:48+00:00,false,N/A,false,N/A,false,N/A
        if not self.password:
            password_enabled = "false"
            password_last_used = "not_supported"
        else:
            password_enabled = "true"
            password_last_used = "no_information"
            if self.password_last_used:
                password_last_used = self.password_last_used.strftime(date_format)

        if len(self.access_keys) == 0:
            access_key_1_active = "false"
            access_key_1_last_rotated = "N/A"
            access_key_1_last_used = "N/A"
            access_key_2_active = "false"
            access_key_2_last_rotated = "N/A"
            access_key_2_last_used = "N/A"
        elif len(self.access_keys) == 1:
            access_key_1_active = (
                "true" if self.access_keys[0].status == "Active" else "false"
            )
            access_key_1_last_rotated = self.access_keys[0].create_date.strftime(
                date_format
            )
            access_key_1_last_used = (
                "N/A"
                if self.access_keys[0].last_used is None
                else self.access_keys[0].last_used.strftime(date_format)
            )
            access_key_2_active = "false"
            access_key_2_last_rotated = "N/A"
            access_key_2_last_used = "N/A"
        else:
            access_key_1_active = (
                "true" if self.access_keys[0].status == "Active" else "false"
            )
            access_key_1_last_rotated = self.access_keys[0].create_date.strftime(
                date_format
            )
            access_key_1_last_used = (
                "N/A"
                if self.access_keys[0].last_used is None
                else self.access_keys[0].last_used.strftime(date_format)
            )
            access_key_2_active = (
                "true" if self.access_keys[1].status == "Active" else "false"
            )
            access_key_2_last_rotated = self.access_keys[1].create_date.strftime(
                date_format
            )
            access_key_2_last_used = (
                "N/A"
                if self.access_keys[1].last_used is None
                else self.access_keys[1].last_used.strftime(date_format)
            )

        cert1_active = cert2_active = False
        if len(self.signing_certificates) > 0:
            cert1 = list(self.signing_certificates.values())[0]
            cert1_active = cert1.status == "Active"
        if len(self.signing_certificates) > 1:
            cert2 = list(self.signing_certificates.values())[1]
            cert2_active = cert2.status == "Active"
        fields = [
            self.name,
            self.arn,
            date_created.strftime(date_format),
            password_enabled,
            password_last_used,
            date_created.strftime(date_format),
            "not_supported",
            "true" if len(self.mfa_devices) else "false",
            access_key_1_active,
            access_key_1_last_rotated,
            access_key_1_last_used,
            "not_supported",
            "not_supported",
            access_key_2_active,
            access_key_2_last_rotated,
            access_key_2_last_used,
            "not_supported",
            "not_supported",
            "true" if cert1_active else "false",
            "N/A",
            "true" if cert2_active else "false",
            "N/A",
        ]
        return ",".join(fields) + "\n"

    @staticmethod
    def cloudformation_name_type() -> str:
        return "UserName"

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::IAM::User"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> User:
        properties = cloudformation_json.get("Properties", {})
        path = properties.get("Path")
        return iam_backends[account_id][get_partition(region_name)].create_user(
            region_name=region_name, user_name=resource_name, path=path
        )

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: User,
        new_resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> User:
        properties = cloudformation_json["Properties"]

        if cls.is_replacement_update(properties):
            resource_name_property = cls.cloudformation_name_type()
            if resource_name_property not in properties:
                properties[resource_name_property] = new_resource_name
            new_resource = cls.create_from_cloudformation_json(
                properties[resource_name_property],
                cloudformation_json,
                account_id,
                region_name,
            )
            properties[resource_name_property] = original_resource.name
            cls.delete_from_cloudformation_json(
                original_resource.name, cloudformation_json, account_id, region_name
            )
            return new_resource

        else:  # No Interruption
            if "Path" in properties:
                original_resource.path = properties["Path"]
            return original_resource

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> None:
        iam_backends[account_id][get_partition(region_name)].delete_user(resource_name)

    @staticmethod
    def is_replacement_update(properties: list[str]) -> bool:
        properties_requiring_replacement_update = ["UserName"]
        return any(
            property_requiring_replacement in properties
            for property_requiring_replacement in properties_requiring_replacement_update
        )

    @property
    def physical_resource_id(self) -> str:
        return self.name


class AccountPasswordPolicy(BaseModel):
    def __init__(
        self,
        allow_change_password: bool,
        hard_expiry: int,
        max_password_age: int,
        minimum_password_length: int,
        password_reuse_prevention: int,
        require_lowercase_characters: bool,
        require_numbers: bool,
        require_symbols: bool,
        require_uppercase_characters: bool,
    ):
        self._errors: list[str] = []
        self._validate(
            max_password_age, minimum_password_length, password_reuse_prevention
        )

        self.allow_users_to_change_password = allow_change_password
        self.hard_expiry = hard_expiry
        self.max_password_age = max_password_age
        self.minimum_password_length = minimum_password_length
        self.password_reuse_prevention = password_reuse_prevention
        self.require_lowercase_characters = require_lowercase_characters
        self.require_numbers = require_numbers
        self.require_symbols = require_symbols
        self.require_uppercase_characters = require_uppercase_characters

    @property
    def expire_passwords(self) -> bool:
        return True if self.max_password_age and self.max_password_age > 0 else False

    def _validate(
        self,
        max_password_age: int,
        minimum_password_length: int,
        password_reuse_prevention: int,
    ) -> None:
        if minimum_password_length > 128:
            self._errors.append(
                self._format_error(
                    key="minimumPasswordLength",
                    value=minimum_password_length,
                    constraint="Member must have value less than or equal to 128",
                )
            )

        if password_reuse_prevention and password_reuse_prevention > 24:
            self._errors.append(
                self._format_error(
                    key="passwordReusePrevention",
                    value=password_reuse_prevention,
                    constraint="Member must have value less than or equal to 24",
                )
            )

        if max_password_age and max_password_age > 1095:
            self._errors.append(
                self._format_error(
                    key="maxPasswordAge",
                    value=max_password_age,
                    constraint="Member must have value less than or equal to 1095",
                )
            )

        self._raise_errors()

    def _format_error(self, key: str, value: Union[str, int], constraint: str) -> str:
        return f'Value "{value}" at "{key}" failed to satisfy constraint: {constraint}'

    def _raise_errors(self) -> None:
        if self._errors:
            count = len(self._errors)
            plural = "s" if len(self._errors) > 1 else ""
            errors = "; ".join(self._errors)
            self._errors = []  # reset collected errors

            raise ValidationError(
                f"{count} validation error{plural} detected: {errors}"
            )


class AccountSummary(BaseModel):
    def __init__(self, iam_backend: IAMBackend):
        self._iam_backend = iam_backend

        self._group_policy_size_quota = 5120
        self._instance_profiles_quota = 1000
        self._groups_per_user_quota = 10
        self._attached_policies_per_user_quota = 10
        self._policies_quota = 1500
        self._account_mfa_enabled = 0  # Haven't found any information being able to activate MFA for the root account programmatically
        self._access_keys_per_user_quota = 2
        self._assume_role_policy_size_quota = 2048
        self._policy_versions_in_use_quota = 10000
        self._global_endpoint_token_version = (
            1  # ToDo: Implement set_security_token_service_preferences()
        )
        self._versions_per_policy_quota = 5
        self._attached_policies_per_group_quota = 10
        self._policy_size_quota = 6144
        self._account_signing_certificates_present = 0  # valid values: 0 | 1
        self._users_quota = 5000
        self._server_certificates_quota = 20
        self._user_policy_size_quota = 2048
        self._roles_quota = 1000
        self._signing_certificates_per_user_quota = 2
        self._role_policy_size_quota = 10240
        self._attached_policies_per_role_quota = 10
        self._account_access_keys_present = 0  # valid values: 0 | 1
        self._groups_quota = 300

    @property
    def summary_map(self) -> dict[str, Any]:  # type: ignore[misc]
        return {
            "GroupPolicySizeQuota": self._group_policy_size_quota,
            "InstanceProfilesQuota": self._instance_profiles_quota,
            "Policies": self._policies,
            "GroupsPerUserQuota": self._groups_per_user_quota,
            "InstanceProfiles": self._instance_profiles,
            "AttachedPoliciesPerUserQuota": self._attached_policies_per_user_quota,
            "Users": self._users,
            "PoliciesQuota": self._policies_quota,
            "Providers": self._providers,
            "AccountMFAEnabled": self._account_mfa_enabled,
            "AccessKeysPerUserQuota": self._access_keys_per_user_quota,
            "AssumeRolePolicySizeQuota": self._assume_role_policy_size_quota,
            "PolicyVersionsInUseQuota": self._policy_versions_in_use_quota,
            "GlobalEndpointTokenVersion": self._global_endpoint_token_version,
            "VersionsPerPolicyQuota": self._versions_per_policy_quota,
            "AttachedPoliciesPerGroupQuota": self._attached_policies_per_group_quota,
            "PolicySizeQuota": self._policy_size_quota,
            "Groups": self._groups,
            "AccountSigningCertificatesPresent": self._account_signing_certificates_present,
            "UsersQuota": self._users_quota,
            "ServerCertificatesQuota": self._server_certificates_quota,
            "MFADevices": self._mfa_devices,
            "UserPolicySizeQuota": self._user_policy_size_quota,
            "PolicyVersionsInUse": self._policy_versions_in_use,
            "ServerCertificates": self._server_certificates,
            "Roles": self._roles,
            "RolesQuota": self._roles_quota,
            "SigningCertificatesPerUserQuota": self._signing_certificates_per_user_quota,
            "MFADevicesInUse": self._mfa_devices_in_use,
            "RolePolicySizeQuota": self._role_policy_size_quota,
            "AttachedPoliciesPerRoleQuota": self._attached_policies_per_role_quota,
            "AccountAccessKeysPresent": self._account_access_keys_present,
            "GroupsQuota": self._groups_quota,
        }

    @property
    def _groups(self) -> int:
        return len(self._iam_backend.groups)

    @property
    def _instance_profiles(self) -> int:
        return len(self._iam_backend.instance_profiles)

    @property
    def _mfa_devices(self) -> int:
        # Don't know, if hardware devices are also counted here
        return len(self._iam_backend.virtual_mfa_devices)

    @property
    def _mfa_devices_in_use(self) -> int:
        devices = 0

        for user in self._iam_backend.users.values():
            devices += len(user.mfa_devices)

        return devices

    @property
    def _policies(self) -> int:
        customer_policies = [
            policy
            for policy in self._iam_backend.managed_policies
            if not re.match(ARN_PARTITION_REGEX + ":iam::aws:policy", policy)
        ]
        return len(customer_policies)

    @property
    def _policy_versions_in_use(self) -> int:
        attachments = 0

        for policy in self._iam_backend.managed_policies.values():
            attachments += policy.attachment_count

        return attachments

    @property
    def _providers(self) -> int:
        return len(self._iam_backend.saml_providers) + len(
            self._iam_backend.open_id_providers
        )

    @property
    def _roles(self) -> int:
        return len(self._iam_backend.roles)

    @property
    def _server_certificates(self) -> int:
        return len(self._iam_backend.certificates)

    @property
    def _users(self) -> int:
        return len(self._iam_backend.users)


def filter_items_with_path_prefix(
    path_prefix: str, items: Iterable[Any]
) -> Iterable[Any]:
    return [role for role in items if role.path.startswith(path_prefix)]


class IAMBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name=region_name, account_id=account_id)
        self.instance_profiles: dict[str, InstanceProfile] = {}
        self.roles: dict[str, Role] = {}
        self.certificates: dict[str, Certificate] = {}
        self.groups: dict[str, Group] = {}
        self.users: dict[str, User] = {}
        self.credential_report: Optional[bool] = None
        self.aws_managed_policies = self._init_aws_policies()
        self.managed_policies = self._init_managed_policies()
        self.account_aliases: list[str] = []
        self.saml_providers: dict[str, SAMLProvider] = {}
        self.open_id_providers: dict[str, OpenIDConnectProvider] = {}
        self.policy_arn_regex = re.compile(
            ARN_PARTITION_REGEX + r":iam::(aws|[0-9]*):policy/.*$"
        )
        self.virtual_mfa_devices: dict[str, VirtualMfaDevice] = {}
        self.account_password_policy: Optional[AccountPasswordPolicy] = None
        self.account_summary = AccountSummary(self)
        self.inline_policies: dict[str, InlinePolicy] = {}
        self.access_keys: dict[str, AccessKey] = {}

        self.tagger = TaggingService()

        self.initialize_service_roles()

    def _init_aws_policies(self) -> list[ManagedPolicy]:
        if not load_iam_aws_managed_policies():
            return []
        # AWS defines some of its own managed policies
        # we periodically import them via `make aws_managed_policies`
        aws_managed_policies_data_parsed = json.loads(aws_managed_policies_data)
        return [
            AWSManagedPolicy.from_data(name, self.account_id, self.region_name, d)
            for name, d in aws_managed_policies_data_parsed.items()
        ]

    def _init_managed_policies(self) -> dict[str, ManagedPolicy]:
        return {p.arn: p for p in self.aws_managed_policies}

    def initialize_service_roles(self) -> None:
        pass
        # TODO: This role is required for some TF tests to work
        # Enabling it breaks an assumption that no roles exist unless created by the user
        #    Our tests, and probably users' tests, rely on this assumption
        # Maybe we can enable this (and roles for other services) as part of a major release
        # self.create_service_linked_role(
        #    service_name="opensearchservice.amazonaws.com", suffix="", description=""
        #    service_name="lakeformation.amazonaws.com"
        # )

    def attach_role_policy(self, policy_arn: str, role_name: str) -> None:
        arns = {p.arn: p for p in self.managed_policies.values()}
        try:
            policy = arns[policy_arn]
        except KeyError:
            raise NotFoundException(
                f"Policy {policy_arn} does not exist or is not attachable."
            )

        policy.attach_to(self.get_role(role_name))

    def update_role_description(self, role_name: str, role_description: str) -> Role:
        role = self.get_role(role_name)
        role.description = role_description
        return role

    def update_role(
        self, role_name: str, role_description: str, max_session_duration: str
    ) -> Role:
        role = self.get_role(role_name)
        role.description = role_description
        role.max_session_duration = max_session_duration
        return role

    def put_role_permissions_boundary(
        self, role_name: str, permissions_boundary: str
    ) -> None:
        if permissions_boundary and not self.policy_arn_regex.match(
            permissions_boundary
        ):
            raise RESTError(
                "InvalidParameterValue",
                f"Value ({permissions_boundary}) for parameter PermissionsBoundary is invalid.",
            )
        role = self.get_role(role_name)
        role.permissions_boundary_arn = permissions_boundary

    def delete_role_permissions_boundary(self, role_name: str) -> None:
        role = self.get_role(role_name)
        role.permissions_boundary_arn = None

    def detach_role_policy(self, policy_arn: str, role_name: str) -> None:
        arns = {p.arn: p for p in self.managed_policies.values()}
        try:
            policy = arns[policy_arn]
            if policy.arn not in self.get_role(role_name).managed_policies.keys():
                raise KeyError
        except KeyError:
            raise NotFoundException(f"Policy {policy_arn} was not found.")
        policy.detach_from(self.get_role(role_name))

    def attach_group_policy(self, policy_arn: str, group_name: str) -> None:
        arns = {p.arn: p for p in self.managed_policies.values()}
        try:
            policy = arns[policy_arn]
        except KeyError:
            raise NotFoundException(f"Policy {policy_arn} was not found.")
        if policy.arn in self.get_group(group_name).managed_policies.keys():
            return
        policy.attach_to(self.get_group(group_name))

    def detach_group_policy(self, policy_arn: str, group_name: str) -> None:
        arns = {p.arn: p for p in self.managed_policies.values()}
        try:
            policy = arns[policy_arn]
            if policy.arn not in self.get_group(group_name).managed_policies.keys():
                raise KeyError
        except KeyError:
            raise NotFoundException(f"Policy {policy_arn} was not found.")
        policy.detach_from(self.get_group(group_name))

    def attach_user_policy(self, policy_arn: str, user_name: str) -> None:
        arns = {p.arn: p for p in self.managed_policies.values()}
        try:
            policy = arns[policy_arn]
        except KeyError:
            raise NotFoundException(
                f"Policy {policy_arn} does not exist or is not attachable."
            )
        policy.attach_to(self.get_user(user_name))

    def detach_user_policy(self, policy_arn: str, user_name: str) -> None:
        arns = {p.arn: p for p in self.managed_policies.values()}
        try:
            policy = arns[policy_arn]
            if policy.arn not in self.get_user(user_name).managed_policies.keys():
                raise KeyError
        except KeyError:
            raise NotFoundException(f"Policy {policy_arn} was not found.")
        policy.detach_from(self.get_user(user_name))

    def create_policy(
        self,
        description: str,
        path: str,
        policy_document: str,
        policy_name: str,
        tags: list[dict[str, str]],
    ) -> ManagedPolicy:
        iam_policy_document_validator = IAMPolicyDocumentValidator(policy_document)
        iam_policy_document_validator.validate()

        clean_tags = self._tag_verification(tags)
        policy = ManagedPolicy(
            policy_name,
            account_id=self.account_id,
            region=self.region_name,
            description=description,
            document=policy_document,
            path=path,
            tags=clean_tags,
        )
        if policy.arn in self.managed_policies:
            raise EntityAlreadyExists(
                f"A policy called {policy_name} already exists. Duplicate names are not allowed."
            )
        self.managed_policies[policy.arn] = policy
        return policy

    def get_policy(self, policy_arn: str) -> ManagedPolicy:
        if policy_arn not in self.managed_policies:
            raise NotFoundException(f"Policy {policy_arn} not found")
        return self.managed_policies[policy_arn]

    def list_attached_role_policies(
        self,
        role_name: str,
        marker: Optional[str] = None,
        max_items: int = 100,
        path_prefix: str = "/",
    ) -> tuple[Iterable[ManagedPolicy], Optional[str]]:
        policies = self.get_role(role_name).managed_policies.values()
        return self._filter_attached_policies(policies, marker, max_items, path_prefix)

    def list_attached_group_policies(
        self,
        group_name: str,
        marker: Optional[str] = None,
        max_items: int = 100,
        path_prefix: str = "/",
    ) -> tuple[Iterable[Any], Optional[str]]:
        policies = self.get_group(group_name).managed_policies.values()
        return self._filter_attached_policies(policies, marker, max_items, path_prefix)

    def list_attached_user_policies(
        self,
        user_name: str,
        marker: Optional[str] = None,
        max_items: int = 100,
        path_prefix: str = "/",
    ) -> tuple[Iterable[Any], Optional[str]]:
        policies = self.get_user(user_name).managed_policies.values()
        return self._filter_attached_policies(policies, marker, max_items, path_prefix)

    def list_policies(
        self,
        marker: Optional[str],
        max_items: int,
        only_attached: bool,
        path_prefix: str,
        scope: str,
    ) -> tuple[Iterable[ManagedPolicy], Optional[str]]:
        policies = list(self.managed_policies.values())

        if only_attached:
            policies = [p for p in policies if p.attachment_count > 0]

        if scope == "AWS":
            policies = [p for p in policies if isinstance(p, AWSManagedPolicy)]
        elif scope == "Local":
            policies = [p for p in policies if not isinstance(p, AWSManagedPolicy)]

        return self._filter_attached_policies(policies, marker, max_items, path_prefix)

    def set_default_policy_version(self, policy_arn: str, version_id: str) -> bool:
        if re.match(r"v[1-9][0-9]*(\.[A-Za-z0-9-]*)?", version_id) is None:
            raise ValidationError(
                f"Value '{version_id}' at 'versionId' failed to satisfy constraint: Member must satisfy regular expression pattern: v[1-9][0-9]*(\\.[A-Za-z0-9-]*)?"
            )

        policy = self.get_policy(policy_arn)

        for version in policy.versions:
            if version.version_id == version_id:
                policy.update_default_version(version_id)
                return True

        raise NoSuchEntity(
            f"Policy {policy_arn} version {version_id} does not exist or is not attachable."
        )

    def _filter_attached_policies(
        self,
        policies: Iterable[Any],
        marker: Optional[str],
        max_items: int,
        path_prefix: str,
    ) -> tuple[Iterable[Any], Optional[str]]:
        if path_prefix:
            policies = [p for p in policies if p.path.startswith(path_prefix)]

        policies = sorted(policies, key=lambda policy: policy.name)
        start_idx = int(marker) if marker else 0

        policies = policies[start_idx : start_idx + max_items]

        if len(policies) < max_items:
            marker = None
        else:
            marker = str(start_idx + max_items)

        return policies, marker

    def create_role(
        self,
        role_name: str,
        assume_role_policy_document: str,
        path: str,
        permissions_boundary: Optional[str],
        description: str,
        tags: list[dict[str, str]],
        max_session_duration: Optional[str],
        linked_service: Optional[str] = None,
    ) -> Role:
        role_id = random_role_id(self.account_id)
        if permissions_boundary and not self.policy_arn_regex.match(
            permissions_boundary
        ):
            raise RESTError(
                "InvalidParameterValue",
                f"Value ({permissions_boundary}) for parameter PermissionsBoundary is invalid.",
            )
        if [role for role in self.get_roles() if role.name == role_name]:
            raise EntityAlreadyExists(f"Role with name {role_name} already exists.")

        clean_tags = self._tag_verification(tags)
        role = Role(
            account_id=self.account_id,
            partition=self.partition,
            role_id=role_id,
            name=role_name,
            assume_role_policy_document=assume_role_policy_document,
            path=path,
            permissions_boundary=permissions_boundary,
            description=description,
            tags=[clean_tags[tag] for tag in clean_tags],
            max_session_duration=max_session_duration,
            linked_service=linked_service,
        )
        self.roles[role_id] = role
        return role

    def get_role_by_id(self, role_id: str) -> Optional[Role]:
        return self.roles.get(role_id)

    def get_role(self, role_name: str) -> Role:
        for role in self.get_roles():
            if role.name == role_name:
                return role
        raise NotFoundException(f"Role {role_name} not found")

    def get_role_by_arn(self, arn: str) -> Role:
        for role in self.get_roles():
            if role.arn == arn:
                return role
        raise NotFoundException(f"Role {arn} not found")

    def delete_role(self, role_name: str) -> None:
        role = self.get_role(role_name)
        for instance_profile in self.get_instance_profiles():
            for profile_role in instance_profile.roles:
                if profile_role.name == role_name:
                    raise DeleteConflictException(
                        "Cannot delete entity, must remove roles from instance profile first."
                    )
        if role.managed_policies:
            raise DeleteConflictException(
                "Cannot delete entity, must detach all policies first."
            )
        if role.policies:
            raise DeleteConflictException(
                "Cannot delete entity, must delete policies first."
            )
        del self.roles[role.id]

    def get_roles(self) -> Iterable[Role]:
        return self.roles.values()

    def update_assume_role_policy(self, role_name: str, policy_document: str) -> None:
        role = self.get_role(role_name)
        iam_policy_document_validator = IAMTrustPolicyDocumentValidator(policy_document)
        iam_policy_document_validator.validate()
        role.assume_role_policy_document = policy_document

    def put_role_policy(
        self, role_name: str, policy_name: str, policy_json: str
    ) -> None:
        role = self.get_role(role_name)

        iam_policy_document_validator = IAMPolicyDocumentValidator(policy_json)
        iam_policy_document_validator.validate()
        role.put_policy(policy_name, policy_json)

    def delete_role_policy(self, role_name: str, policy_name: str) -> None:
        role = self.get_role(role_name)
        role.delete_policy(policy_name)

    def get_role_policy(self, role_name: str, policy_name: str) -> tuple[str, str]:
        role = self.get_role(role_name)
        for p, d in role.policies.items():
            if p == policy_name:
                return p, d
        raise NotFoundException(
            f"Policy Document {policy_name} not attached to role {role_name}"
        )

    def list_role_policies(self, role_name: str) -> list[str]:
        role = self.get_role(role_name)
        return list(role.policies.keys())

    def _tag_verification(
        self, tags: list[dict[str, str]]
    ) -> dict[str, dict[str, str]]:
        if len(tags) > 50:
            raise TooManyTags(tags)

        tag_keys: dict[str, dict[str, str]] = {}
        for tag in tags:
            # Need to index by the lowercase tag key since the keys are case insensitive, but their case is retained.
            ref_key = tag["Key"].lower()
            self._check_tag_duplicate(tag_keys, ref_key)
            self._validate_tag_key(tag["Key"])
            if len(tag["Value"]) > 256:
                raise TagValueTooBig(tag["Value"])

            tag_keys[ref_key] = tag

        return tag_keys

    def _validate_tag_key(
        self, tag_key: str, exception_param: str = "tags.X.member.key"
    ) -> None:
        """Validates the tag key.

        :param tag_key: The tag key to check against.
        :param exception_param: The exception parameter to send over to help format the message. This is to reflect
                                the difference between the tag and untag APIs.
        :return:
        """
        # Validate that the key length is correct:
        if len(tag_key) > 128:
            raise TagKeyTooBig(tag_key, param=exception_param)

        # Validate that the tag key fits the proper Regex:
        # [\w\s_.:/=+\-@]+ SHOULD be the same as the Java regex on the AWS documentation: [\p{L}\p{Z}\p{N}_.:/=+\-@]+
        match = re.findall(r"[\w\s_.:/=+\-@]+", tag_key)
        # Kudos if you can come up with a better way of doing a global search :)
        if not len(match) or len(match[0]) < len(tag_key):
            raise InvalidTagCharacters(tag_key, param=exception_param)

    def _check_tag_duplicate(
        self, all_tags: dict[str, dict[str, str]], tag_key: str
    ) -> None:
        """Validates that a tag key is not a duplicate

        :param all_tags: Dict to check if there is a duplicate tag.
        :param tag_key: The tag key to check against.
        :return:
        """
        if tag_key in all_tags:
            raise DuplicateTags()

    def list_role_tags(
        self, role_name: str, marker: Optional[str], max_items: int = 100
    ) -> tuple[list[dict[str, str]], Optional[str]]:
        role = self.get_role(role_name)

        max_items = int(max_items)
        tag_index = sorted(role.tags, key=lambda t: t["Key"])
        start_idx = int(marker) if marker else 0

        tag_index = tag_index[start_idx : start_idx + max_items]

        if len(role.tags) <= (start_idx + max_items):
            marker = None
        else:
            marker = str(start_idx + max_items)

        return tag_index, marker

    def tag_role(self, role_name: str, tags: list[dict[str, str]]) -> None:
        self._tag_verification(tags)
        role = self.get_role(role_name)
        new_keys = [tag["Key"] for tag in tags]
        updated_tags = [tag for tag in role.tags if tag["Key"] not in new_keys]
        updated_tags.extend(tags)
        role.tags = updated_tags

    def untag_role(self, role_name: str, tag_keys: list[str]) -> None:
        if len(tag_keys) > 50:
            raise TooManyTags(tag_keys, param="tagKeys")

        role = self.get_role(role_name)

        for key in tag_keys:
            self._validate_tag_key(key, exception_param="tagKeys")

        role.tags = [tag for tag in role.tags if tag["Key"] not in tag_keys]

    def list_policy_tags(
        self, policy_arn: str, marker: Optional[str], max_items: int = 100
    ) -> tuple[list[dict[str, str]], Optional[str]]:
        policy = self.get_policy(policy_arn)

        max_items = int(max_items)
        tag_index = sorted(policy.tags)
        start_idx = int(marker) if marker else 0

        tag_index = tag_index[start_idx : start_idx + max_items]

        if len(policy.tags) <= (start_idx + max_items):
            marker = None
        else:
            marker = str(start_idx + max_items)

        # Make the tag list of dict's:
        tags = [policy.tags[tag] for tag in tag_index]

        return tags, marker

    def tag_policy(self, policy_arn: str, tags: list[dict[str, str]]) -> None:
        clean_tags = self._tag_verification(tags)
        policy = self.get_policy(policy_arn)
        policy.tags.update(clean_tags)

    def untag_policy(self, policy_arn: str, tag_keys: list[str]) -> None:
        if len(tag_keys) > 50:
            raise TooManyTags(tag_keys, param="tagKeys")

        policy = self.get_policy(policy_arn)

        for key in tag_keys:
            ref_key = key.lower()
            self._validate_tag_key(key, exception_param="tagKeys")

            policy.tags.pop(ref_key, None)

    def create_policy_version(
        self, policy_arn: str, policy_document: str, set_as_default: bool
    ) -> PolicyVersion:
        iam_policy_document_validator = IAMPolicyDocumentValidator(policy_document)
        iam_policy_document_validator.validate()

        policy = self.get_policy(policy_arn)
        if not policy:
            raise NotFoundException("Policy not found")
        if len(policy.versions) >= 5:
            raise LimitExceededException(
                "A managed policy can have up to 5 versions. Before you create a new version, you must delete an existing version."
            )
        version = PolicyVersion(policy_arn, policy_document, set_as_default)
        policy.versions.append(version)
        version.version_id = f"v{policy.next_version_num}"
        policy.next_version_num += 1
        if set_as_default:
            policy.update_default_version(version.version_id)
        return version

    def get_policy_version(self, policy_arn: str, version_id: str) -> PolicyVersion:
        policy = self.get_policy(policy_arn)
        if not policy:
            raise NotFoundException("Policy not found")
        for version in policy.versions:
            if version.version_id == version_id:
                return version
        raise NotFoundException("Policy version not found")

    def list_policy_versions(self, policy_arn: str) -> list[PolicyVersion]:
        policy = self.get_policy(policy_arn)
        if not policy:
            raise NotFoundException("Policy not found")
        return policy.versions

    def delete_policy_version(self, policy_arn: str, version_id: str) -> None:
        policy = self.get_policy(policy_arn)
        if not policy:
            raise NotFoundException("Policy not found")
        if version_id == policy.default_version_id:
            raise DeleteConflictException(
                "Cannot delete the default version of a policy."
            )
        for i, v in enumerate(policy.versions):
            if v.version_id == version_id:
                del policy.versions[i]
                return
        raise NotFoundException("Policy not found")

    def create_instance_profile(
        self,
        name: str,
        path: str,
        role_names: list[str],
        tags: Optional[list[dict[str, str]]] = None,
    ) -> InstanceProfile:
        if self.instance_profiles.get(name):
            raise EntityAlreadyExists(
                f"Instance Profile {name} already exists.",
            )

        instance_profile_id = random_resource_id()

        roles = [self.get_role(role_name) for role_name in role_names]
        instance_profile = InstanceProfile(
            account_id=self.account_id,
            region_name=self.region_name,
            instance_profile_id=instance_profile_id,
            name=name,
            path=path,
            roles=roles,
            tags=tags,
        )
        self.instance_profiles[name] = instance_profile
        return instance_profile

    def delete_instance_profile(
        self, name: str, ignore_attached_roles: bool = False
    ) -> None:
        instance_profile = self.get_instance_profile(name)
        if len(instance_profile.roles) > 0 and not ignore_attached_roles:
            raise DeleteConflictException(
                "Cannot delete entity, must remove roles from instance profile first.",
            )
        del self.instance_profiles[name]

    def get_instance_profile(self, profile_name: str) -> InstanceProfile:
        for profile in self.get_instance_profiles():
            if profile.name == profile_name:
                return profile

        raise NotFoundException(f"Instance profile {profile_name} not found")

    def get_instance_profile_by_arn(self, profile_arn: str) -> InstanceProfile:
        for profile in self.get_instance_profiles():
            if profile.arn == profile_arn:
                return profile

        raise NotFoundException(f"Instance profile {profile_arn} not found")

    def get_instance_profiles(self) -> Iterable[InstanceProfile]:
        return self.instance_profiles.values()

    def get_instance_profiles_for_role(self, role_name: str) -> list[InstanceProfile]:
        found_profiles = []

        for profile in self.get_instance_profiles():
            if len(profile.roles) > 0:
                if profile.roles[0].name == role_name:
                    found_profiles.append(profile)

        return found_profiles

    def add_role_to_instance_profile(self, profile_name: str, role_name: str) -> None:
        profile = self.get_instance_profile(profile_name)
        role = self.get_role(role_name)
        if not profile.roles:
            profile.roles.append(role)
        else:
            raise LimitExceededException(
                "Cannot exceed quota for InstanceSessionsPerInstanceProfile: 1"
            )

    def remove_role_from_instance_profile(
        self, profile_name: str, role_name: str
    ) -> None:
        profile = self.get_instance_profile(profile_name)
        role = self.get_role(role_name)
        profile.roles.remove(role)

    def list_server_certificates(self) -> Iterable[Certificate]:
        """
        Pagination is not yet implemented
        """
        return self.certificates.values()

    def upload_server_certificate(
        self,
        cert_name: str,
        cert_body: str,
        private_key: str,
        cert_chain: Optional[str] = None,
        path: Optional[str] = None,
    ) -> Certificate:
        certificate_id = random_resource_id()
        cert = Certificate(
            account_id=self.account_id,
            region_name=self.region_name,
            cert_name=cert_name,
            cert_body=cert_body,
            private_key=private_key,
            cert_chain=cert_chain,
            path=path,
        )
        self.certificates[certificate_id] = cert
        return cert

    def get_server_certificate(self, name: str) -> Certificate:
        for cert in self.certificates.values():
            if name == cert.cert_name:
                return cert

        raise NotFoundException(
            f"The Server Certificate with name {name} cannot be found."
        )

    def get_certificate_by_arn(self, arn: str) -> Optional[Certificate]:
        for cert in self.certificates.values():
            if arn == cert.arn:
                return cert
        return None

    def delete_server_certificate(self, name: str) -> None:
        cert_id = None
        for key, cert in self.certificates.items():
            if name == cert.cert_name:
                cert_id = key
                break

        if cert_id is None:
            raise NotFoundException(
                f"The Server Certificate with name {name} cannot be found."
            )

        self.certificates.pop(cert_id, None)

    def create_group(self, group_name: str, path: str = "/") -> Group:
        if group_name in self.groups:
            raise EntityAlreadyExists(f"Group {group_name} already exists")

        group = Group(self.account_id, self.region_name, group_name, path)
        self.groups[group_name] = group
        return group

    def get_group(self, group_name: str) -> Group:
        """
        Pagination is not yet implemented
        """
        try:
            return self.groups[group_name]
        except KeyError:
            raise NotFoundException(f"Group {group_name} not found")

    def list_groups(self) -> Iterable[Group]:
        return self.groups.values()

    def get_groups_for_user(self, user_name: str) -> list[Group]:
        user = self.get_user(user_name)
        groups = []
        for group in self.list_groups():
            if user in group.users:
                groups.append(group)

        return groups

    def put_group_policy(
        self, group_name: str, policy_name: str, policy_json: str
    ) -> None:
        group = self.get_group(group_name)

        iam_policy_document_validator = IAMPolicyDocumentValidator(policy_json)
        iam_policy_document_validator.validate()
        group.put_policy(policy_name, policy_json)

    def list_group_policies(self, group_name: str) -> list[str]:
        """
        Pagination is not yet implemented
        """
        group = self.get_group(group_name)
        return group.list_policies()

    def delete_group_policy(self, group_name: str, policy_name: str) -> None:
        group = self.get_group(group_name)
        group.delete_policy(policy_name)

    def get_group_policy(self, group_name: str, policy_name: str) -> dict[str, str]:
        group = self.get_group(group_name)
        return group.get_policy(policy_name)

    def delete_group(self, group_name: str) -> None:
        try:
            del self.groups[group_name]
        except KeyError:
            raise NotFoundException(
                f"The group with name {group_name} cannot be found."
            )

    def update_group(
        self, group_name: str, new_group_name: Optional[str], new_path: Optional[str]
    ) -> None:
        if new_group_name:
            if new_group_name in self.groups:
                raise EntityAlreadyExists(f"Group {new_group_name} already exists")
            try:
                group = self.groups[group_name]
            except KeyError:
                raise NotFoundException(
                    f"The group with name {group_name} cannot be found."
                )

            existing_policies = group.managed_policies.copy()
            for policy_arn in existing_policies:
                self.detach_group_policy(policy_arn, group_name)
            if new_path:
                group.path = new_path
            group.name = new_group_name
            self.groups[new_group_name] = self.groups.pop(group_name)
            for policy_arn in existing_policies:
                self.attach_group_policy(policy_arn, new_group_name)

    def create_user(
        self,
        region_name: str,
        user_name: str,
        path: str = "/",
        tags: Optional[list[dict[str, str]]] = None,
    ) -> User:
        if user_name in self.users:
            raise EntityAlreadyExists(f"User {user_name} already exists")

        user = User(self.account_id, region_name, user_name, path)
        self.tagger.tag_resource(user.arn, tags or [])
        self.users[user_name] = user
        return user

    def get_user(self, name: str) -> User:
        user = self.users.get(name)

        if not user:
            raise NoSuchEntity(f"The user with name {name} cannot be found.")

        return user

    def list_users(
        self,
        path_prefix: Optional[str],
        marker: Optional[str],
        max_items: Optional[int],
    ) -> Iterable[User]:
        try:
            users: Iterable[User] = list(self.users.values())
            if path_prefix:
                users = filter_items_with_path_prefix(path_prefix, users)

        except KeyError:
            raise NotFoundException(
                f"Users {path_prefix}, {marker}, {max_items} not found"
            )

        return users

    def update_user(
        self,
        user_name: str,
        new_path: Optional[str] = None,
        new_user_name: Optional[str] = None,
    ) -> None:
        try:
            user = self.users[user_name]
        except KeyError:
            raise NotFoundException(f"User {user_name} not found")

        if new_path:
            user.path = new_path
        if new_user_name:
            user.name = new_user_name
            self.users[new_user_name] = self.users.pop(user_name)

    def list_roles(
        self,
        path_prefix: Optional[str] = None,
        marker: Optional[str] = None,
        max_items: Optional[int] = None,
    ) -> tuple[list[Role], Optional[str]]:
        path_prefix = path_prefix if path_prefix else "/"
        max_items = int(max_items) if max_items else 100
        start_index = int(marker) if marker else 0

        roles: Iterable[Role] = list(self.roles.values())
        roles = filter_items_with_path_prefix(path_prefix, roles)
        sorted_roles = sorted(roles, key=lambda role: role.id)

        roles_to_return = sorted_roles[start_index : start_index + max_items]

        if len(sorted_roles) <= (start_index + max_items):
            marker = None
        else:
            marker = str(start_index + max_items)

        return roles_to_return, marker

    def upload_signing_certificate(
        self, user_name: str, body: str
    ) -> SigningCertificate:
        user = self.get_user(user_name)
        cert_id = random_resource_id(size=32)

        # Validate the signing cert:
        try:
            data = bytes(body, "utf8")

            x509.load_pem_x509_certificate(data, default_backend())

        except Exception:
            raise MalformedCertificate(body)

        if (
            len(user.signing_certificates)
            >= self.account_summary._signing_certificates_per_user_quota
        ):
            raise LimitExceededException(
                "Cannot exceed quota for CertificatesPerUser: 2"
            )
        user.signing_certificates[cert_id] = SigningCertificate(
            cert_id, user_name, body
        )

        return user.signing_certificates[cert_id]

    def delete_signing_certificate(self, user_name: str, cert_id: str) -> None:
        user = self.get_user(user_name)

        try:
            del user.signing_certificates[cert_id]
        except KeyError:
            raise NotFoundException(
                f"The Certificate with id {cert_id} cannot be found."
            )

    def list_signing_certificates(self, user_name: str) -> list[SigningCertificate]:
        user = self.get_user(user_name)

        return list(user.signing_certificates.values())

    def update_signing_certificate(
        self, user_name: str, cert_id: str, status: str
    ) -> None:
        user = self.get_user(user_name)

        try:
            user.signing_certificates[cert_id].status = status

        except KeyError:
            raise NotFoundException(
                f"The Certificate with id {cert_id} cannot be found."
            )

    def create_login_profile(self, user_name: str, password: str) -> User:
        # This does not currently deal with PasswordPolicyViolation.
        user = self.get_user(user_name)
        if user.password:
            raise EntityAlreadyExists(f"User {user_name} already has password")
        user.password = password
        return user

    def get_login_profile(self, user_name: str) -> User:
        user = self.get_user(user_name)
        if not user.password:
            raise NotFoundException(f"Login profile for {user_name} not found")
        return user

    def update_login_profile(
        self, user_name: str, password: str, password_reset_required: bool
    ) -> User:
        # This does not currently deal with PasswordPolicyViolation.
        user = self.get_user(user_name)
        if not user.password:
            raise NotFoundException(f"Login profile for {user_name} not found")
        user.password = password
        user.password_reset_required = password_reset_required
        return user

    def delete_login_profile(self, user_name: str) -> None:
        user = self.get_user(user_name)
        if not user.password:
            raise NotFoundException(f"Login profile for {user_name} not found")
        user.password = None

    def add_user_to_group(self, group_name: str, user_name: str) -> None:
        user = self.get_user(user_name)
        group = self.get_group(group_name)
        if user not in group.users:
            group.users.append(user)

    def remove_user_from_group(self, group_name: str, user_name: str) -> None:
        group = self.get_group(group_name)
        user = self.get_user(user_name)
        try:
            group.users.remove(user)
        except ValueError:
            raise NotFoundException(f"User {user_name} not in group {group_name}")

    def get_user_policy(self, user_name: str, policy_name: str) -> dict[str, str]:
        user = self.get_user(user_name)
        return user.get_policy(policy_name)

    def list_user_policies(self, user_name: str) -> Iterable[str]:
        user = self.get_user(user_name)
        return user.policies.keys()

    def list_user_tags(self, user_name: str) -> list[dict[str, str]]:
        user = self.get_user(user_name)
        return self.tagger.list_tags_for_resource(user.arn)["Tags"]

    def put_user_policy(
        self, user_name: str, policy_name: str, policy_json: str
    ) -> None:
        user = self.get_user(user_name)

        iam_policy_document_validator = IAMPolicyDocumentValidator(policy_json)
        iam_policy_document_validator.validate()
        user.put_policy(policy_name, policy_json)

    def delete_user_policy(self, user_name: str, policy_name: str) -> None:
        user = self.get_user(user_name)
        user.delete_policy(policy_name)

    def delete_policy(self, policy_arn: str) -> None:
        policy = self.get_policy(policy_arn)
        del self.managed_policies[policy.arn]

    def create_access_key(
        self, user_name: str, prefix: str = "AKIA", status: str = "Active"
    ) -> AccessKey:
        keys = self.list_access_keys(user_name)
        if len(keys) >= LIMIT_KEYS_PER_USER:
            raise LimitExceededException(
                f"Cannot exceed quota for AccessKeysPerUser: {LIMIT_KEYS_PER_USER}"
            )
        user = self.get_user(user_name)
        key = user.create_access_key(prefix=prefix, status=status)
        self.access_keys[key.physical_resource_id] = key
        return key

    def create_temp_access_key(self) -> AccessKey:
        # Temporary access keys such as the ones returned by STS when assuming a role temporarily
        key = AccessKey(user_name=None, prefix="ASIA", account_id=self.account_id)

        self.access_keys[key.physical_resource_id] = key
        return key

    def update_access_key(
        self, user_name: str, access_key_id: str, status: Optional[str] = None
    ) -> AccessKey:
        user = self.get_user(user_name)
        return user.update_access_key(access_key_id, status)

    def get_access_key_last_used(self, access_key_id: str) -> dict[str, Any]:
        access_keys_list = self.get_all_access_keys_for_all_users()
        for key in access_keys_list:
            if key.access_key_id == access_key_id:
                return {"user_name": key.user_name, "last_used": key.last_used}

        raise NotFoundException(
            f"The Access Key with id {access_key_id} cannot be found"
        )

    def get_all_access_keys_for_all_users(self) -> list[AccessKey]:
        access_keys_list = []
        for account in iam_backends.values():
            for user_name in account[self.partition].users:
                access_keys_list += account[self.partition].list_access_keys(user_name)
        return access_keys_list

    def list_access_keys(self, user_name: str) -> list[AccessKey]:
        """
        Pagination is not yet implemented
        """
        user = self.get_user(user_name)
        return user.get_all_access_keys()

    def delete_access_key(self, access_key_id: str, user_name: str) -> None:
        user = self.get_user(user_name)
        access_key = user.get_access_key_by_id(access_key_id)
        self.delete_access_key_by_name(access_key.access_key_id)

    def delete_access_key_by_name(self, name: str) -> None:
        key = self.access_keys[name]
        try:  # User may have been deleted before their access key...
            if key.user_name is not None:
                user = self.get_user(key.user_name)
                user.delete_access_key(key.access_key_id)
        except NoSuchEntity:
            pass
        del self.access_keys[name]

    def upload_ssh_public_key(
        self, user_name: str, ssh_public_key_body: str
    ) -> SshPublicKey:
        user = self.get_user(user_name)
        return user.upload_ssh_public_key(ssh_public_key_body)

    def get_ssh_public_key(
        self, user_name: str, ssh_public_key_id: str
    ) -> SshPublicKey:
        user = self.get_user(user_name)
        return user.get_ssh_public_key(ssh_public_key_id)

    def get_all_ssh_public_keys(self, user_name: str) -> Iterable[SshPublicKey]:
        user = self.get_user(user_name)
        return user.get_all_ssh_public_keys()

    def update_ssh_public_key(
        self, user_name: str, ssh_public_key_id: str, status: str
    ) -> None:
        user = self.get_user(user_name)
        user.update_ssh_public_key(ssh_public_key_id, status)

    def delete_ssh_public_key(self, user_name: str, ssh_public_key_id: str) -> None:
        user = self.get_user(user_name)
        user.delete_ssh_public_key(ssh_public_key_id)

    def enable_mfa_device(
        self,
        user_name: str,
        serial_number: str,
        authentication_code_1: str,
        authentication_code_2: str,
    ) -> None:
        """Enable MFA Device for user."""
        user = self.get_user(user_name)
        if serial_number in user.mfa_devices:
            raise EntityAlreadyExists(f"Device {serial_number} already exists")

        device = self.virtual_mfa_devices.get(serial_number, None)
        if device:
            device.enable_date = utcnow()
            device.user = user
            device.user_attribute = {
                "Path": user.path,
                "UserName": user.name,
                "UserId": user.id,
                "Arn": user.arn,
                "CreateDate": user.create_date,
                "PasswordLastUsed": None,  # not supported
                "PermissionsBoundary": {},  # ToDo: add put_user_permissions_boundary() functionality
                "Tags": self.tagger.list_tags_for_resource(user.arn)["Tags"],
            }

        user.enable_mfa_device(
            serial_number, authentication_code_1, authentication_code_2
        )

    def deactivate_mfa_device(self, user_name: str, serial_number: str) -> None:
        """Deactivate and detach MFA Device from user if device exists."""
        user = self.get_user(user_name)
        if serial_number not in user.mfa_devices:
            raise NotFoundException(f"Device {serial_number} not found")

        device = self.virtual_mfa_devices.get(serial_number, None)
        if device:
            device.enable_date = None
            device.user = None
            device.user_attribute = None

        user.deactivate_mfa_device(serial_number)

    def list_mfa_devices(self, user_name: str) -> Iterable[MFADevice]:
        user = self.get_user(user_name)
        return user.mfa_devices.values()

    def create_virtual_mfa_device(
        self, device_name: str, path: str
    ) -> VirtualMfaDevice:
        if not path:
            path = "/"

        if not path.startswith("/") and not path.endswith("/"):
            raise ValidationError(
                "The specified value for path is invalid. "
                "It must begin and end with / and contain only alphanumeric characters and/or / characters."
            )

        if any(not len(part) for part in path.split("/")[1:-1]):
            raise ValidationError(
                "The specified value for path is invalid. "
                "It must begin and end with / and contain only alphanumeric characters and/or / characters."
            )

        if len(path) > 512:
            raise ValidationError(
                "1 validation error detected: "
                'Value "{}" at "path" failed to satisfy constraint: '
                "Member must have length less than or equal to 512"
            )

        device = VirtualMfaDevice(
            self.account_id,
            region_name=self.region_name,
            device_name=path + device_name,
        )

        if device.serial_number in self.virtual_mfa_devices:
            raise EntityAlreadyExists(
                "MFADevice entity at the same path and name already exists."
            )

        self.virtual_mfa_devices[device.serial_number] = device
        return device

    def delete_virtual_mfa_device(self, serial_number: str) -> None:
        device = self.virtual_mfa_devices.pop(serial_number, None)

        if not device:
            raise NotFoundException(
                f"VirtualMFADevice with serial number {serial_number} doesn't exist."
            )

    def list_virtual_mfa_devices(
        self, assignment_status: str, marker: Optional[str], max_items: int
    ) -> tuple[list[VirtualMfaDevice], Optional[str]]:
        devices = list(self.virtual_mfa_devices.values())

        if assignment_status == "Assigned":
            devices = [device for device in devices if device.enable_date]

        if assignment_status == "Unassigned":
            devices = [device for device in devices if not device.enable_date]

        sorted(devices, key=lambda device: device.serial_number)
        max_items = int(max_items)
        start_idx = int(marker) if marker else 0

        if start_idx > len(devices):
            raise ValidationError("Invalid Marker.")

        devices = devices[start_idx : start_idx + max_items]

        if len(devices) < max_items:
            marker = None
        else:
            marker = str(start_idx + max_items)

        return devices, marker

    def delete_user(self, user_name: str) -> None:
        user = self.get_user(user_name)
        if user.managed_policies:
            raise DeleteConflictException(
                "Cannot delete entity, must detach all policies first."
            )
        if user.policies:
            raise DeleteConflictException(
                "Cannot delete entity, must delete policies first."
            )
        self.tagger.delete_all_tags_for_resource(user.arn)
        del self.users[user_name]

    def report_generated(self) -> Optional[bool]:
        return self.credential_report

    def generate_report(self) -> None:
        self.credential_report = True

    def get_credential_report(self) -> str:
        if not self.credential_report:
            raise ReportNotPresentException("Credential report not present")
        report = "user,arn,user_creation_time,password_enabled,password_last_used,password_last_changed,password_next_rotation,mfa_active,access_key_1_active,access_key_1_last_rotated,access_key_1_last_used_date,access_key_1_last_used_region,access_key_1_last_used_service,access_key_2_active,access_key_2_last_rotated,access_key_2_last_used_date,access_key_2_last_used_region,access_key_2_last_used_service,cert_1_active,cert_1_last_rotated,cert_2_active,cert_2_last_rotated\n"
        for user in self.users:
            report += self.users[user].to_csv()
        return report

    def list_account_aliases(self) -> list[str]:
        return self.account_aliases

    def create_account_alias(self, alias: str) -> None:
        # alias is force updated
        self.account_aliases = [alias]

    def delete_account_alias(self) -> None:
        self.account_aliases = []

    def get_account_authorization_details(
        self, policy_filter: list[str]
    ) -> dict[str, Any]:
        policies = self.managed_policies.values()
        local_policies = set(policies) - set(self.aws_managed_policies)
        returned_policies = []

        if len(policy_filter) == 0:
            return {
                "instance_profiles": self.instance_profiles.values(),
                "roles": self.roles.values(),
                "groups": self.groups.values(),
                "users": self.users.values(),
                "managed_policies": self.managed_policies.values(),
            }

        if "AWSManagedPolicy" in policy_filter:
            returned_policies = self.aws_managed_policies
        if "LocalManagedPolicy" in policy_filter:
            returned_policies = returned_policies + list(local_policies)

        return {
            "instance_profiles": self.instance_profiles.values(),
            "roles": self.roles.values() if "Role" in policy_filter else [],
            "groups": self.groups.values() if "Group" in policy_filter else [],
            "users": self.users.values() if "User" in policy_filter else [],
            "managed_policies": returned_policies,
        }

    def create_saml_provider(
        self, name: str, saml_metadata_document: str
    ) -> SAMLProvider:
        saml_provider = SAMLProvider(
            account_id=self.account_id,
            region_name=self.region_name,
            name=name,
            saml_metadata_document=saml_metadata_document,
        )
        self.saml_providers[name] = saml_provider
        return saml_provider

    def update_saml_provider(
        self, saml_provider_arn: str, saml_metadata_document: str
    ) -> SAMLProvider:
        saml_provider = self.get_saml_provider(saml_provider_arn)
        saml_provider.saml_metadata_document = saml_metadata_document
        return saml_provider

    def delete_saml_provider(self, saml_provider_arn: str) -> None:
        try:
            for saml_provider in list(self.list_saml_providers()):
                if saml_provider.arn == saml_provider_arn:
                    del self.saml_providers[saml_provider.name]
        except KeyError:
            raise NotFoundException(f"SAMLProvider {saml_provider_arn} not found")

    def list_saml_providers(self) -> Iterable[SAMLProvider]:
        return self.saml_providers.values()

    def get_saml_provider(self, saml_provider_arn: str) -> SAMLProvider:
        for saml_provider in self.list_saml_providers():
            if saml_provider.arn == saml_provider_arn:
                return saml_provider
        raise NotFoundException(f"SamlProvider {saml_provider_arn} not found")

    def get_user_from_access_key_id(self, access_key_id: str) -> Optional[User]:
        for user_name, user in self.users.items():
            access_keys = self.list_access_keys(user_name)
            for access_key in access_keys:
                if access_key.access_key_id == access_key_id:
                    return user
        return None

    def create_open_id_connect_provider(
        self,
        url: str,
        thumbprint_list: list[str],
        client_id_list: list[str],
        tags: list[dict[str, str]],
    ) -> OpenIDConnectProvider:
        clean_tags = self._tag_verification(tags)
        open_id_provider = OpenIDConnectProvider(
            account_id=self.account_id,
            region_name=self.region_name,
            url=url,
            thumbprint_list=thumbprint_list,
            client_id_list=client_id_list,
            tags=clean_tags,
        )

        if open_id_provider.arn in self.open_id_providers:
            raise EntityAlreadyExists("Unknown")

        self.open_id_providers[open_id_provider.arn] = open_id_provider
        return open_id_provider

    def update_open_id_connect_provider_thumbprint(
        self, arn: str, thumbprint_list: list[str]
    ) -> None:
        open_id_provider = self.get_open_id_connect_provider(arn)
        open_id_provider.thumbprint_list = thumbprint_list

    def tag_open_id_connect_provider(
        self, arn: str, tags: list[dict[str, str]]
    ) -> None:
        open_id_provider = self.get_open_id_connect_provider(arn)
        clean_tags = self._tag_verification(tags)
        open_id_provider.tags.update(clean_tags)

    def untag_open_id_connect_provider(self, arn: str, tag_keys: list[str]) -> None:
        open_id_provider = self.get_open_id_connect_provider(arn)

        for key in tag_keys:
            ref_key = key.lower()
            self._validate_tag_key(key, exception_param="tagKeys")
            open_id_provider.tags.pop(ref_key, None)

    def list_open_id_connect_provider_tags(
        self, arn: str, marker: Optional[str], max_items: int = 100
    ) -> tuple[list[dict[str, str]], Optional[str]]:
        open_id_provider = self.get_open_id_connect_provider(arn)

        max_items = int(max_items)
        tag_index = sorted(open_id_provider.tags)
        start_idx = int(marker) if marker else 0

        tag_index = tag_index[start_idx : start_idx + max_items]

        if len(open_id_provider.tags) <= (start_idx + max_items):
            marker = None
        else:
            marker = str(start_idx + max_items)

        tags = [open_id_provider.tags[tag] for tag in tag_index]
        return tags, marker

    def delete_open_id_connect_provider(self, arn: str) -> None:
        self.open_id_providers.pop(arn, None)

    def get_open_id_connect_provider(self, arn: str) -> OpenIDConnectProvider:
        open_id_provider = self.open_id_providers.get(arn)

        if not open_id_provider:
            raise NotFoundException(f"OpenIDConnect Provider not found for arn {arn}")

        return open_id_provider

    def list_open_id_connect_providers(self) -> list[str]:
        return list(self.open_id_providers.keys())

    def update_account_password_policy(
        self,
        allow_change_password: bool,
        hard_expiry: int,
        max_password_age: int,
        minimum_password_length: int,
        password_reuse_prevention: int,
        require_lowercase_characters: bool,
        require_numbers: bool,
        require_symbols: bool,
        require_uppercase_characters: bool,
    ) -> None:
        self.account_password_policy = AccountPasswordPolicy(
            allow_change_password,
            hard_expiry,
            max_password_age,
            minimum_password_length,
            password_reuse_prevention,
            require_lowercase_characters,
            require_numbers,
            require_symbols,
            require_uppercase_characters,
        )

    def get_account_password_policy(self) -> AccountPasswordPolicy:
        if not self.account_password_policy:
            raise NoSuchEntity(
                f"The Password Policy with domain name {self.account_id} cannot be found."
            )

        return self.account_password_policy

    def delete_account_password_policy(self) -> None:
        if not self.account_password_policy:
            raise NoSuchEntity(
                "The account policy with name PasswordPolicy cannot be found."
            )

        self.account_password_policy = None

    def get_account_summary(self) -> AccountSummary:
        return self.account_summary

    def create_inline_policy(
        self,
        resource_name: str,
        policy_name: str,
        policy_document: str,
        group_names: list[str],
        role_names: list[str],
        user_names: list[str],
    ) -> InlinePolicy:
        if resource_name in self.inline_policies:
            raise EntityAlreadyExists(f"Inline Policy {resource_name} already exists")

        inline_policy = InlinePolicy(
            resource_name,
            policy_name,
            policy_document,
            group_names,
            role_names,
            user_names,
        )
        self.inline_policies[resource_name] = inline_policy
        inline_policy.apply_policy(self)
        return inline_policy

    def get_inline_policy(self, policy_id: str) -> InlinePolicy:
        try:
            return self.inline_policies[policy_id]
        except KeyError:
            raise NotFoundException(f"Inline policy {policy_id} not found")

    def update_inline_policy(
        self,
        resource_name: str,
        policy_name: str,
        policy_document: str,
        group_names: list[str],
        role_names: list[str],
        user_names: list[str],
    ) -> InlinePolicy:
        inline_policy = self.get_inline_policy(resource_name)
        inline_policy.unapply_policy(self)
        inline_policy.update(
            policy_name, policy_document, group_names, role_names, user_names
        )
        inline_policy.apply_policy(self)
        return inline_policy

    def delete_inline_policy(self, policy_id: str) -> None:
        inline_policy = self.get_inline_policy(policy_id)
        inline_policy.unapply_policy(self)
        del self.inline_policies[policy_id]

    def tag_user(self, name: str, tags: list[dict[str, str]]) -> None:
        user = self.get_user(name)

        self.tagger.tag_resource(user.arn, tags)

    def untag_user(self, name: str, tag_keys: list[str]) -> None:
        user = self.get_user(name)

        self.tagger.untag_resource_using_names(user.arn, tag_keys)

    def create_service_linked_role(
        self, service_name: str, description: str, suffix: str
    ) -> Role:
        # service.amazonaws.com -> Service
        # some-thing.service.amazonaws.com -> Service_SomeThing
        service = service_name.split(".")[-3]
        prefix = service_name.split(".")[0]
        if service != prefix:
            prefix = "".join([x.capitalize() for x in prefix.split("-")])
            service = SERVICE_NAME_CONVERSION.get(service, service) + "_" + prefix
        else:
            service = SERVICE_NAME_CONVERSION.get(service, service)
        role_name = f"AWSServiceRoleFor{service}"
        if suffix:
            role_name = role_name + f"_{suffix}"
        assume_role_policy_document = {
            "Version": "2012-10-17",
            "Statement": [
                {
                    "Action": ["sts:AssumeRole"],
                    "Effect": "Allow",
                    "Principal": {"Service": [service_name]},
                }
            ],
        }
        path = f"/aws-service-role/{service_name}/"
        return self.create_role(
            role_name,
            json.dumps(assume_role_policy_document),
            path,
            permissions_boundary=None,
            description=description,
            tags=[],
            max_session_duration="3600",
            linked_service=service_name,
        )

    def delete_service_linked_role(self, role_name: str) -> str:
        self.delete_role(role_name)
        deletion_task_id = str(random.uuid4())
        return deletion_task_id

    def get_service_linked_role_deletion_status(self) -> bool:
        """
        This method always succeeds for now - we do not yet keep track of deletions
        """
        return True

    def tag_instance_profile(
        self, instance_profile_name: str, tags: list[dict[str, str]]
    ) -> None:
        profile = self.get_instance_profile(profile_name=instance_profile_name)
        new_keys = [tag["Key"] for tag in tags]
        updated_tags = [tag for tag in profile.tags if tag["Key"] not in new_keys]
        updated_tags.extend(tags)
        profile.tags = updated_tags

    def untag_instance_profile(
        self, instance_profile_name: str, tagKeys: list[str]
    ) -> None:
        profile = self.get_instance_profile(profile_name=instance_profile_name)
        profile.tags = [tag for tag in profile.tags if tag["Key"] not in tagKeys]


iam_backends = BackendDict(
    IAMBackend, "iam", use_boto3_regions=False, additional_regions=PARTITION_NAMES
)
