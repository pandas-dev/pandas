from typing import Any, Dict, Iterable

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random as random
from moto.utilities.utils import get_partition

from .exceptions import ResourceNotFoundException


def _create_id(aws_account_id: str, namespace: str, _id: str) -> str:
    return f"{aws_account_id}:{namespace}:{_id}"


class QuicksightDataSet(BaseModel):
    def __init__(self, account_id: str, region: str, _id: str, name: str):
        self.arn = f"arn:{get_partition(region)}:quicksight:{region}:{account_id}:data-set/{_id}"
        self._id = _id
        self.name = name
        self.region = region
        self.account_id = account_id

    def to_json(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "DataSetId": self._id,
            "IngestionArn": f"arn:{get_partition(self.region)}:quicksight:{self.region}:{self.account_id}:ingestion/tbd",
        }


class QuicksightIngestion(BaseModel):
    def __init__(
        self, account_id: str, region: str, data_set_id: str, ingestion_id: str
    ):
        self.arn = f"arn:{get_partition(region)}:quicksight:{region}:{account_id}:data-set/{data_set_id}/ingestions/{ingestion_id}"
        self.ingestion_id = ingestion_id

    def to_json(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "IngestionId": self.ingestion_id,
            "IngestionStatus": "INITIALIZED",
        }


class QuicksightMembership(BaseModel):
    def __init__(self, account_id: str, region: str, group: str, user: str):
        self.group = group
        self.user = user
        self.arn = f"arn:{get_partition(region)}:quicksight:{region}:{account_id}:group/default/{group}/{user}"

    def to_json(self) -> Dict[str, str]:
        return {"Arn": self.arn, "MemberName": self.user}


class QuicksightGroup(BaseModel):
    def __init__(
        self,
        region: str,
        group_name: str,
        description: str,
        aws_account_id: str,
        namespace: str,
    ):
        self.arn = f"arn:{get_partition(region)}:quicksight:{region}:{aws_account_id}:group/default/{group_name}"
        self.group_name = group_name
        self.description = description
        self.aws_account_id = aws_account_id
        self.namespace = namespace
        self.region = region

        self.members: Dict[str, QuicksightMembership] = dict()

    def add_member(self, user_name: str) -> QuicksightMembership:
        membership = QuicksightMembership(
            self.aws_account_id, self.region, self.group_name, user_name
        )
        self.members[user_name] = membership
        return membership

    def delete_member(self, user_name: str) -> None:
        self.members.pop(user_name, None)

    def get_member(self, user_name: str) -> QuicksightMembership:
        return self.members[user_name]

    def list_members(self) -> Iterable[QuicksightMembership]:
        return self.members.values()

    def to_json(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "GroupName": self.group_name,
            "Description": self.description,
            "PrincipalId": self.aws_account_id,
            "Namespace": self.namespace,
        }


class QuicksightUser(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        email: str,
        identity_type: str,
        username: str,
        user_role: str,
    ):
        self.arn = f"arn:{get_partition(region)}:quicksight:{region}:{account_id}:user/default/{username}"
        self.email = email
        self.identity_type = identity_type
        self.username = username
        self.user_role = user_role
        self.active = False
        self.principal_id = random.get_random_hex(10)

    def to_json(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "Email": self.email,
            "IdentityType": self.identity_type,
            "Role": self.user_role,
            "UserName": self.username,
            "Active": self.active,
            "PrincipalId": self.principal_id,
        }


class QuickSightBackend(BaseBackend):
    """Implementation of QuickSight APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.groups: Dict[str, QuicksightGroup] = dict()
        self.users: Dict[str, QuicksightUser] = dict()

    def create_data_set(self, data_set_id: str, name: str) -> QuicksightDataSet:
        return QuicksightDataSet(
            self.account_id, self.region_name, data_set_id, name=name
        )

    def create_group(
        self, group_name: str, description: str, aws_account_id: str, namespace: str
    ) -> QuicksightGroup:
        group = QuicksightGroup(
            region=self.region_name,
            group_name=group_name,
            description=description,
            aws_account_id=aws_account_id,
            namespace=namespace,
        )
        _id = _create_id(aws_account_id, namespace, group_name)
        self.groups[_id] = group
        return group

    def create_group_membership(
        self, aws_account_id: str, namespace: str, group_name: str, user_name: str
    ) -> QuicksightMembership:
        group = self.describe_group(aws_account_id, namespace, group_name)
        return group.add_member(user_name)

    def create_ingestion(
        self, data_set_id: str, ingestion_id: str
    ) -> QuicksightIngestion:
        return QuicksightIngestion(
            self.account_id, self.region_name, data_set_id, ingestion_id
        )

    def delete_group(
        self, aws_account_id: str, namespace: str, group_name: str
    ) -> None:
        _id = _create_id(aws_account_id, namespace, group_name)
        self.groups.pop(_id, None)

    def delete_user(self, aws_account_id: str, namespace: str, user_name: str) -> None:
        # Delete users from all groups
        for group in self.groups.values():
            group.delete_member(user_name)
        # Delete user itself
        _id = _create_id(aws_account_id, namespace, user_name)
        self.users.pop(_id, None)

    def describe_group(
        self, aws_account_id: str, namespace: str, group_name: str
    ) -> QuicksightGroup:
        _id = _create_id(aws_account_id, namespace, group_name)
        if _id not in self.groups:
            raise ResourceNotFoundException(f"Group {group_name} not found")
        return self.groups[_id]

    def describe_group_membership(
        self, aws_account_id: str, namespace: str, group_name: str, user_name: str
    ) -> QuicksightMembership:
        group = self.describe_group(aws_account_id, namespace, group_name)
        return group.get_member(user_name)

    def describe_user(
        self, aws_account_id: str, namespace: str, user_name: str
    ) -> QuicksightUser:
        _id = _create_id(aws_account_id, namespace, user_name)
        if _id not in self.users:
            raise ResourceNotFoundException(f"User {user_name} not found")
        return self.users[_id]

    def list_groups(
        self, aws_account_id: str, namespace: str
    ) -> Iterable[QuicksightGroup]:
        """
        The NextToken and MaxResults parameters are not yet implemented
        """
        id_for_ns = _create_id(aws_account_id, namespace, _id="")
        return [
            group for _id, group in self.groups.items() if _id.startswith(id_for_ns)
        ]

    def list_group_memberships(
        self, aws_account_id: str, namespace: str, group_name: str
    ) -> Iterable[QuicksightMembership]:
        """
        The NextToken and MaxResults parameters are not yet implemented
        """
        group = self.describe_group(aws_account_id, namespace, group_name)
        return group.list_members()

    def list_users(
        self, aws_account_id: str, namespace: str
    ) -> Iterable[QuicksightUser]:
        """
        The NextToken and MaxResults parameters are not yet implemented
        """
        id_for_ns = _create_id(aws_account_id, namespace, _id="")
        return [user for _id, user in self.users.items() if _id.startswith(id_for_ns)]

    def register_user(
        self,
        identity_type: str,
        email: str,
        user_role: str,
        aws_account_id: str,
        namespace: str,
        user_name: str,
    ) -> QuicksightUser:
        """
        The following parameters are not yet implemented:
        IamArn, SessionName, CustomsPermissionsName, ExternalLoginFederationProviderType, CustomFederationProviderUrl, ExternalLoginId
        """
        user = QuicksightUser(
            account_id=self.account_id,
            region=self.region_name,
            email=email,
            identity_type=identity_type,
            user_role=user_role,
            username=user_name,
        )
        _id = _create_id(aws_account_id, namespace, user_name)
        self.users[_id] = user
        return user

    def update_group(
        self, aws_account_id: str, namespace: str, group_name: str, description: str
    ) -> QuicksightGroup:
        group = self.describe_group(aws_account_id, namespace, group_name)
        group.description = description
        return group


quicksight_backends = BackendDict(QuickSightBackend, "quicksight")
