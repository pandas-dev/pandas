import re
import string
from typing import Any, Dict, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time, utcnow
from moto.moto_api._internal import mock_random as random
from moto.organizations.models import OrganizationsBackend, organizations_backends
from moto.ram.exceptions import (
    InvalidParameterException,
    MalformedArnException,
    OperationNotPermittedException,
    UnknownResourceException,
)
from moto.utilities.utils import get_partition


def random_resource_id(size: int) -> str:
    return "".join(random.choice(string.digits + "abcdef") for _ in range(size))


class ResourceShare(BaseModel):
    # List of shareable resources can be found here
    # https://docs.aws.amazon.com/ram/latest/userguide/shareable.html
    SHAREABLE_RESOURCES = [
        "cluster",  # Amazon Aurora cluster
        "component",  # Amazon EC2 Image Builder component
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
        "transit-gateway",  # Amazon EC2 transit gateway
    ]

    def __init__(self, account_id: str, region: str, **kwargs: Any):
        self.account_id = account_id
        self.region = region
        self.partition = get_partition(region)

        self.allow_external_principals = kwargs.get("allowExternalPrincipals", True)
        self.arn = f"arn:{get_partition(self.region)}:ram:{self.region}:{account_id}:resource-share/{random.uuid4()}"
        self.creation_time = utcnow()
        self.feature_set = "STANDARD"
        self.last_updated_time = utcnow()
        self.name = kwargs["name"]
        self.owning_account_id = account_id
        self.principals: List[str] = []
        self.resource_arns: List[str] = []
        self.status = "ACTIVE"

    @property
    def organizations_backend(self) -> OrganizationsBackend:
        return organizations_backends[self.account_id][self.partition]

    def add_principals(self, principals: List[str]) -> None:
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
                    if any(principal == ou["Arn"] for ou in ous):
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

    def add_resources(self, resource_arns: List[str]) -> None:
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

    def describe(self) -> Dict[str, Any]:
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

    def update(self, **kwargs: Any) -> None:
        self.allow_external_principals = kwargs.get(
            "allowExternalPrincipals", self.allow_external_principals
        )
        self.last_updated_time = utcnow()
        self.name = kwargs.get("name", self.name)


class ResourceAccessManagerBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.resource_shares: List[ResourceShare] = []

    @property
    def organizations_backend(self) -> OrganizationsBackend:
        return organizations_backends[self.account_id][self.partition]

    def create_resource_share(self, **kwargs: Any) -> Dict[str, Any]:
        resource = ResourceShare(self.account_id, self.region_name, **kwargs)
        resource.add_principals(kwargs.get("principals", []))
        resource.add_resources(kwargs.get("resourceArns", []))

        self.resource_shares.append(resource)

        response = resource.describe()
        response.pop("featureSet")

        return dict(resourceShare=response)

    def get_resource_shares(self, **kwargs: Any) -> Dict[str, Any]:
        owner = kwargs["resourceOwner"]

        if owner not in ["SELF", "OTHER-ACCOUNTS"]:
            raise InvalidParameterException(
                f"{owner} is not a valid resource owner. "
                "Specify either SELF or OTHER-ACCOUNTS and try again."
            )

        if owner == "OTHER-ACCOUNTS":
            raise NotImplementedError(
                "Value 'OTHER-ACCOUNTS' for parameter 'resourceOwner' not implemented."
            )

        resouces = [resource.describe() for resource in self.resource_shares]

        return dict(resourceShares=resouces)

    def update_resource_share(self, **kwargs: Any) -> Dict[str, Any]:
        arn = kwargs["resourceShareArn"]

        resource = next(
            (resource for resource in self.resource_shares if arn == resource.arn), None
        )

        if not resource:
            raise UnknownResourceException(f"ResourceShare {arn} could not be found.")

        resource.update(**kwargs)
        response = resource.describe()
        response.pop("featureSet")

        return dict(resourceShare=response)

    def delete_resource_share(self, arn: str) -> Dict[str, Any]:
        resource = next(
            (resource for resource in self.resource_shares if arn == resource.arn), None
        )

        if not resource:
            raise UnknownResourceException(f"ResourceShare {arn} could not be found.")

        resource.delete()

        return dict(returnValue=True)

    def enable_sharing_with_aws_organization(self) -> Dict[str, Any]:
        if not self.organizations_backend.org:
            raise OperationNotPermittedException

        return dict(returnValue=True)


ram_backends = BackendDict(ResourceAccessManagerBackend, "ram")
