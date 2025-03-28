from typing import Any, Dict, List, Union

from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random as random
from moto.utilities.utils import get_partition


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

    def add_member(self, member_name: str) -> QuicksightMembership:
        membership = QuicksightMembership(
            self.aws_account_id, self.region, self.group_name, member_name
        )
        self.members[member_name] = membership
        return membership

    def delete_member(self, user_name: str) -> None:
        self.members.pop(user_name, None)

    def get_member(self, user_name: str) -> Union[QuicksightMembership, None]:
        return self.members.get(user_name, None)

    def list_members(self) -> List[QuicksightMembership]:
        return list(self.members.values())

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
