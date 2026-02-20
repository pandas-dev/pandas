import re
import string
from datetime import datetime, timezone
from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time, utcnow
from moto.moto_api._internal import mock_random as random
from moto.organizations.models import (
    OrganizationsBackend,
    organizations_backends,
)
from moto.ram.exceptions import (
    InvalidParameterException,
    MalformedArnException,
    OperationNotPermittedException,
    UnknownResourceException,
)
from moto.ram.utils import AWS_MANAGED_PERMISSIONS, RAM_RESOURCE_TYPES
from moto.utilities.utils import get_partition


def random_resource_id(size: int) -> str:
    return "".join(random.choice(string.digits + "abcdef") for _ in range(size))


class ResourceShare(BaseModel):
    # List of shareable resources can be found here
    # https://docs.aws.amazon.com/ram/latest/userguide/shareable.html
    SHAREABLE_RESOURCES = [
        "cluster",  # Amazon Aurora cluster
        "component",  # Amazon EC2 Image Builder component
        "core-network",  # Amazon Network Manager core network
        "group",  # AWS Resource Groups
        "image",  # Amazon EC2 Image Builder image
        "image-recipe",  # Amazon EC2 Image Builder image recipe
        "license-configuration",  # AWS License Manager configuration
        "mesh",  # AWS App Mesh
        "prefix-list",  # Amazon EC2 prefix list
        "project",  # AWS CodeBuild project
        "report-group",  # AWS CodeBuild report group
        "resolver-rule",  # Amazon Route 53 forwarding rule
        "subnet",  # Amazon EC2 subnet
        "transit-gateway",  # Amazon EC2 transit gateway,
        "database",  # Amazon Glue database
        "table",  # Amazon Glue table
        "catalog",  # Amazon Glue catalog
    ]

    def __init__(self, account_id: str, region: str, **kwargs: Any):
        self.account_id = account_id
        self.region = region
        self.partition = get_partition(region)

        self.allow_external_principals = kwargs.get("allowExternalPrincipals", True)
        self.arn = f"arn:{self.partition}:ram:{self.region}:{account_id}:resource-share/{random.uuid4()}"
        self.creation_time = utcnow()
        self.feature_set = "STANDARD"
        self.last_updated_time = utcnow()
        self.name = kwargs["name"]
        self.owning_account_id = account_id
        self.principals: list[str] = []
        self.resource_arns: list[str] = []
        self.status = "ACTIVE"

    @property
    def organizations_backend(self) -> OrganizationsBackend:
        return organizations_backends[self.account_id][self.partition]

    def add_principals(self, principals: list[str]) -> None:
        for principal in principals:
            match = re.search(
                r"^arn:aws:organizations::\d{12}:organization/(o-\w+)$", principal
            )
            if match:
                organization = self.organizations_backend.describe_organization()
                if principal == organization["Organization"]["Arn"]:
                    continue
                else:
                    raise UnknownResourceException(
                        f"Organization {match.group(1)} could not be found."
                    )

            match = re.search(
                r"^arn:aws:organizations::\d{12}:ou/(o-\w+)/(ou-[\w-]+)$", principal
            )
            if match:
                roots = self.organizations_backend.list_roots()
                root_id = next(
                    (
                        root["Id"]
                        for root in roots["Roots"]
                        if root["Name"] == "Root" and match.group(1) in root["Arn"]
                    ),
                    None,
                )

                if root_id:
                    (
                        ous,
                        _,
                    ) = self.organizations_backend.list_organizational_units_for_parent(
                        parent_id=root_id
                    )
                    if any(principal == ou.arn for ou in ous):
                        continue

                raise UnknownResourceException(
                    f"OrganizationalUnit {match.group(2)} in unknown organization could not be found."
                )

            if not re.match(r"^\d{12}$", principal):
                raise InvalidParameterException(
                    f"Principal ID {principal} is malformed. Verify the ID and try again."
                )

        for principal in principals:
            self.principals.append(principal)

    def add_resources(self, resource_arns: list[str]) -> None:
        for resource in resource_arns:
            match = re.search(
                r"^arn:aws:[a-z0-9-]+:[a-z0-9-]*:[0-9]{12}:([a-z-]+)[/:].*$", resource
            )
            if not match:
                raise MalformedArnException(
                    f"The specified resource ARN {resource} is not valid. Verify the ARN and try again."
                )

            if match.group(1) not in self.SHAREABLE_RESOURCES:
                raise MalformedArnException(
                    "You cannot share the selected resource type."
                )

        for resource in resource_arns:
            self.resource_arns.append(resource)

    def delete(self) -> None:
        self.last_updated_time = utcnow()
        self.status = "DELETED"

    def describe(self) -> dict[str, Any]:
        return {
            "allowExternalPrincipals": self.allow_external_principals,
            "creationTime": unix_time(self.creation_time),
            "featureSet": self.feature_set,
            "lastUpdatedTime": unix_time(self.last_updated_time),
            "name": self.name,
            "owningAccountId": self.owning_account_id,
            "resourceShareArn": self.arn,
            "status": self.status,
        }

    def update(self, allow_external_principals: bool, name: Optional[str]) -> None:
        self.allow_external_principals = allow_external_principals
        self.last_updated_time = utcnow()
        self.name = name


class ResourceType(BaseModel):
    def __init__(
        self, resource_type: str, service_name: str, resource_region_scope: str
    ):
        self.resource_type = resource_type
        self.service_name = service_name
        self.resource_region_scope = resource_region_scope

    def describe(self) -> dict[str, Any]:
        return {
            "resourceType": self.resource_type,
            "serviceName": self.service_name,
            "resourceRegionScope": self.resource_region_scope,
        }


class ManagedPermission(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        name: str,
        resource_type: str,
        version: str = "1",
        default_version: bool = True,
        status: str = "ATTACHABLE",
        creation_time: Optional[str] = None,
        last_updated_time: Optional[str] = None,
        is_resource_type_default: bool = False,
        permission_type: str = "AWS_MANAGED",  # or "CUSTOMER_MANAGED",
    ):
        self.account_id = account_id
        self.region = region
        self.partition = get_partition(region)

        now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
        arn_prefix = (
            f"arn:{self.partition}:ram::{self.partition}:permission/"
            if permission_type == "AWS_MANAGED"
            else f"arn:{self.partition}:ram:{self.region}:{account_id}:permission/"
        )

        self.name = name
        self.arn = f"{arn_prefix}{name}"
        self.resource_type = resource_type
        self.version = version
        self.default_version = default_version
        self.status = status
        self.creation_time = creation_time or now
        self.last_updated_time = last_updated_time or creation_time
        self.is_resource_type_default = is_resource_type_default
        self.permission_type = permission_type

    def describe(self) -> dict[str, Any]:
        return {
            "arn": self.arn,
            "name": self.name,
            "resourceType": self.resource_type,
            "version": self.version,
            "defaultVersion": self.default_version,
            "status": self.status,
            "creationTime": self.creation_time,
            "lastUpdatedTime": self.last_updated_time,
            "isResourceTypeDefault": self.is_resource_type_default,
            "permissionType": self.permission_type,
        }


class ResourceAccessManagerBackend(BaseBackend):
    PERMISSION_TYPES = ["ALL", "AWS", "LOCAL"]

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.resource_shares: list[ResourceShare] = []
        self.managed_permissions: list[ManagedPermission] = [
            ManagedPermission(
                account_id=account_id,
                region=region_name,
                name=permission["name"],
                resource_type=permission["resourceType"],
                version=permission["version"],
                default_version=permission["defaultVersion"],
                status=permission["status"],
                creation_time=permission["creationTime"],
                last_updated_time=permission["lastUpdatedTime"],
                is_resource_type_default=permission["isResourceTypeDefault"],
                permission_type=permission["permissionType"],
            )
            for permission in AWS_MANAGED_PERMISSIONS
        ]
        self.resource_types: list[ResourceType] = [
            ResourceType(
                resource_type=resource_type["resourceType"],
                service_name=resource_type["serviceName"],
                resource_region_scope=resource_type["resourceRegionScope"],
            )
            for resource_type in RAM_RESOURCE_TYPES
        ]

    @property
    def organizations_backend(self) -> OrganizationsBackend:
        return organizations_backends[self.account_id][self.partition]

    def create_resource_share(self, **kwargs: Any) -> dict[str, Any]:
        resource = ResourceShare(self.account_id, self.region_name, **kwargs)
        resource.add_principals(kwargs.get("principals", []))
        resource.add_resources(kwargs.get("resourceArns", []))

        self.resource_shares.append(resource)

        response = resource.describe()
        response.pop("featureSet")

        return {"resourceShare": response}

    def get_resource_shares(self, resource_owner: Optional[str]) -> dict[str, Any]:
        if resource_owner not in ["SELF", "OTHER-ACCOUNTS"]:
            raise InvalidParameterException(
                f"{resource_owner} is not a valid resource owner. "
                "Specify either SELF or OTHER-ACCOUNTS and try again."
            )

        if resource_owner == "OTHER-ACCOUNTS":
            raise NotImplementedError(
                "Value 'OTHER-ACCOUNTS' for parameter 'resourceOwner' not implemented."
            )

        resources = [resource.describe() for resource in self.resource_shares]

        return {"resourceShares": resources}

    def update_resource_share(
        self,
        resource_share_arn: Optional[str],
        allow_external_principals: bool,
        name: Optional[str],
    ) -> dict[str, Any]:
        resource = next(
            (
                resource
                for resource in self.resource_shares
                if resource_share_arn == resource.arn
            ),
            None,
        )

        if not resource:
            raise UnknownResourceException(
                f"ResourceShare {resource_share_arn} could not be found."
            )

        resource.update(allow_external_principals, name)
        response = resource.describe()
        response.pop("featureSet")

        return {"resourceShare": response}

    def delete_resource_share(
        self, resource_share_arn: Optional[str]
    ) -> dict[str, Any]:
        resource = next(
            (
                resource
                for resource in self.resource_shares
                if resource_share_arn == resource.arn
            ),
            None,
        )

        if not resource:
            raise UnknownResourceException(
                f"ResourceShare {resource_share_arn} could not be found."
            )

        resource.delete()

        return {"returnValue": True}

    def enable_sharing_with_aws_organization(self) -> dict[str, Any]:
        if not self.organizations_backend.org:
            raise OperationNotPermittedException

        return {"returnValue": True}

    def get_resource_share_associations(
        self,
        association_type: Optional[str],
        association_status: Optional[str],
        resource_share_arns: list[str],
        resource_arn: Optional[str],
        principal: Optional[str],
    ) -> dict[str, Any]:
        if association_type not in ["PRINCIPAL", "RESOURCE"]:
            raise InvalidParameterException(
                f"{association_type} is not a valid association type. "
                "Specify either PRINCIPAL or RESOURCE and try again."
            )

        if association_status and association_status not in [
            "ASSOCIATING",
            "ASSOCIATED",
            "FAILED",
            "DISASSOCIATING",
            "DISASSOCIATED",
        ]:
            raise InvalidParameterException(
                f"{association_status} is not a valid association status."
            )

        if association_type == "PRINCIPAL" and resource_arn:
            raise InvalidParameterException(
                "You cannot specify a resource ARN when the association type is PRINCIPAL."
            )
        if association_type == "RESOURCE" and principal:
            raise InvalidParameterException(
                "You cannot specify a principal when the association type is RESOURCE."
            )

        associations = []
        for resource_share in self.resource_shares:
            if resource_share_arns and resource_share.arn not in resource_share_arns:
                continue

            if association_type == "PRINCIPAL":
                for principal_id in resource_share.principals:
                    if principal and principal != principal_id:
                        continue
                    associations.append(
                        {
                            "resourceShareArn": resource_share.arn,
                            "resourceShareName": resource_share.name,
                            "associatedEntity": principal_id,
                            "associationType": "PRINCIPAL",
                            "status": association_status or "ASSOCIATED",
                            "creationTime": unix_time(resource_share.creation_time),
                            "lastUpdatedTime": unix_time(
                                resource_share.last_updated_time
                            ),
                            "external": False,
                        }
                    )
            else:  # RESOURCE
                for resource_id in resource_share.resource_arns:
                    if resource_arn and resource_arn != resource_id:
                        continue
                    associations.append(
                        {
                            "resourceShareArn": resource_share.arn,
                            "resourceShareName": resource_share.name,
                            "associatedEntity": resource_id,
                            "associationType": "RESOURCE",
                            "status": association_status or "ASSOCIATED",
                            "creationTime": unix_time(resource_share.creation_time),
                            "lastUpdatedTime": unix_time(
                                resource_share.last_updated_time
                            ),
                            "external": False,
                        }
                    )

        return {"resourceShareAssociations": associations}

    def list_resource_types(self, resource_region_scope: str) -> dict[str, Any]:
        if resource_region_scope not in ["ALL", "REGIONAL", "GLOBAL"]:
            raise InvalidParameterException(
                f"{resource_region_scope} is not a valid resource region "
                "scope value. Specify a valid value and try again."
            )

        if resource_region_scope == "ALL":
            resource_types = [
                resource_type.describe() for resource_type in self.resource_types
            ]
        else:
            resource_types = [
                resource_type_dict
                for resource_type in self.resource_types
                if (resource_type_dict := resource_type.describe())
                and resource_type_dict["resourceRegionScope"] == resource_region_scope
            ]

        return {"resourceTypes": resource_types}

    def list_permissions(
        self, resource_type: str, permission_type: str
    ) -> dict[str, Any]:
        permission_types_relation = {
            "AWS": "AWS_MANAGED",
            "LOCAL": "CUSTOMER_MANAGED",
        }

        # Here, resourceType first partition (service) is case sensitive and
        # last partition (type) is case insensitive
        if resource_type and not any(
            (permission_dict := permission.describe())
            and permission_dict["resourceType"].split(":")[0]
            == resource_type.split(":")[0]
            and permission_dict["resourceType"].lower() == resource_type.lower()
            for permission in self.managed_permissions
        ):
            raise InvalidParameterException(f"Invalid resource type: {resource_type}")

        if resource_type:
            permissions = [
                permission_dict
                for permission in self.managed_permissions
                if (permission_dict := permission.describe())
                and permission_dict["resourceType"].lower() == resource_type.lower()
            ]
        else:
            permissions = [
                permission.describe() for permission in self.managed_permissions
            ]

        if permission_type not in self.PERMISSION_TYPES:
            raise InvalidParameterException(
                f"{permission_type} is not a valid scope. Must be one of: "
                f"{', '.join(self.PERMISSION_TYPES)}."
            )

        if permission_type != "ALL":
            permissions = [
                permission
                for permission in permissions
                if permission_types_relation.get(permission_type)
                == permission["permissionType"]
            ]

        return {"permissions": permissions}


ram_backends = BackendDict(ResourceAccessManagerBackend, "ram")
