import datetime
from typing import Any, Optional, Union

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

    def to_json(self) -> dict[str, Any]:
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

    def to_json(self) -> dict[str, Any]:
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

    def to_json(self) -> dict[str, str]:
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

        self.members: dict[str, QuicksightMembership] = {}

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

    def list_members(self) -> list[QuicksightMembership]:
        return list(self.members.values())

    def to_json(self) -> dict[str, Any]:
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

    def to_json(self) -> dict[str, Any]:
        return {
            "Arn": self.arn,
            "Email": self.email,
            "IdentityType": self.identity_type,
            "Role": self.user_role,
            "UserName": self.username,
            "Active": self.active,
            "PrincipalId": self.principal_id,
        }


class QuicksightDashboard(BaseModel):
    # Matches model from https://docs.aws.amazon.com/quicksight/latest/APIReference/API_CreateDashboard.html#API_CreateDashboard_RequestSyntax
    # Todo: Handle versions and all fields
    def __init__(
        self,
        account_id: str,
        region: str,
        dashboard_id: str,
        dashboard_publish_options: dict[str, Any],
        name: str,
        definition: dict[str, Any],
        folder_arns: list[str],
        link_entities: list[str],
        link_sharing_configuration: dict[str, Any],
        parameters: dict[str, Any],
        permissions: list[dict[str, Any]],
        source_entity: dict[str, Any],
        tags: list[dict[str, Any]],
        theme_arn: str,
        version_description: str,
        validation_strategy: dict[str, str],
    ) -> None:
        self.arn = f"arn:{get_partition(region)}:quicksight:{region}:{account_id}:dashboard/{dashboard_id}"
        self.dashboard_id = dashboard_id
        self.name = name
        self.region = region
        self.account_id = account_id
        self.dashboard_publish_options = dashboard_publish_options
        self.definition = definition
        self.folder_arns = folder_arns
        self.link_entities = link_entities
        self.link_sharing_configuration = link_sharing_configuration
        self.parameters = parameters
        self.permissions = permissions
        self.source_entity = source_entity
        self.tags = tags
        self.theme_arn = theme_arn
        self.version_description = version_description
        self.validation_strategy = validation_strategy
        # Not user provided
        self.created_time = datetime.datetime.now()
        self.version_number = 1
        self.status = "CREATION_SUCCESSFUL"
        self.last_updated_time = datetime.datetime.now()
        self.last_published_time = datetime.datetime.now()

    def to_dict(self) -> dict[str, Any]:
        return {
            "Arn": self.arn,
            "DashboardId": self.dashboard_id,
            "Name": self.name,
            "Version": {
                "CreatedTime": str(self.created_time),
                "Errors": [],
                "Status": self.status,
                "VersionNumber": self.version_number,
                "Arn": self.arn,
                "SourceEntityArn": self.source_entity,
                "ThemeArn": self.theme_arn,
                "Description": self.version_description,
                "SourceEntity": self.source_entity,
            },
            "CreatedTime": str(self.created_time),
            "LastPublishedTime": str(self.last_published_time),
            "LastUpdatedTime": str(self.last_updated_time),
        }


class QuicksightAccountSettings(BaseModel):
    def __init__(
        self, account_id: str, account_name: Optional[str] = "default"
    ) -> None:
        self.account_name = account_name
        self.account_id = account_id
        self.default_namespace = "default"
        self.notification_email = ""
        self.termination_protection_enabled = False
        self.public_sharing_enabled = False
        self.edition = "STANDARD"


class QuickSightDataSource(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        data_source_id: str,
        name: str,
        data_source_parameters: Optional[dict[str, dict[str, Any]]] = None,
        alternate_data_source_parameters: Optional[list[dict[str, Any]]] = None,
        ssl_properties: Optional[dict[str, Any]] = None,
        status: Optional[str] = None,
        tags: Optional[list[dict[str, str]]] = None,
        data_source_type: Optional[str] = None,
        vpc_connection_properties: Optional[dict[str, Any]] = None,
    ) -> None:
        self.account_id = account_id
        self.region = region
        self.created_time = datetime.datetime.now()
        self.data_source_id = data_source_id
        self.data_source_parameters = data_source_parameters
        self.alternate_data_source_parameters = alternate_data_source_parameters
        self.last_updated_time = datetime.datetime.now()
        self.name = name
        self.data_source_type = data_source_type
        self.arn = f"arn:{get_partition(region)}:quicksight:{region}:{account_id}:datasource/{data_source_id}"
        self.status = status
        self.tags = tags or []
        self.ssl_properties = ssl_properties
        self.vpc_connection_properties = vpc_connection_properties

    def to_json(self) -> dict[str, Any]:
        return {
            "AlternateDataSourceParameters": self.alternate_data_source_parameters,
            "Arn": self.arn,
            "CreatedTime": self.created_time.isoformat(),
            "DataSourceId": self.data_source_id,
            "DataSourceParameters": self.data_source_parameters,
            "LastUpdatedTime": self.last_updated_time.isoformat(),
            "Name": self.name,
            "SslProperties": self.ssl_properties,
            "Status": self.status,
            "Type": self.data_source_type,
            "VpcConnectionProperties": self.vpc_connection_properties,
        }
