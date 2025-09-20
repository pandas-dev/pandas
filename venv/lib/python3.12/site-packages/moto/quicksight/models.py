from datetime import datetime
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService

from .data_models import (
    QuicksightAccountSettings,
    QuicksightDashboard,
    QuicksightDataSet,
    QuickSightDataSource,
    QuicksightGroup,
    QuicksightIngestion,
    QuicksightMembership,
    QuicksightUser,
)
from .exceptions import ResourceNotFoundException
from .utils import PAGINATION_MODEL, QuicksightSearchFilterFactory


def _create_id(aws_account_id: str, namespace: str, _id: str) -> str:
    return f"{aws_account_id}:{namespace}:{_id}"


class QuickSightBackend(BaseBackend):
    """Implementation of QuickSight APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.dashboards: Dict[str, QuicksightDashboard] = dict()
        self.groups: Dict[str, QuicksightGroup] = dict()
        self.users: Dict[str, QuicksightUser] = dict()
        self.account_settings: QuicksightAccountSettings = QuicksightAccountSettings(
            account_id=account_id
        )
        self.data_sources: Dict[str, QuickSightDataSource] = dict()
        self.data_sets: Dict[str, QuicksightDataSet] = dict()
        self.tagger = TaggingService()

    def create_data_set(
        self,
        data_set_id: str,
        name: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> QuicksightDataSet:
        dataset = QuicksightDataSet(
            self.account_id, self.region_name, data_set_id, name=name
        )

        if tags:
            self.tagger.tag_resource(
                arn=dataset.arn,
                tags=tags,
            )

        self.data_sets[data_set_id] = dataset
        return dataset

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
        self, aws_account_id: str, namespace: str, group_name: str, member_name: str
    ) -> QuicksightMembership:
        group = self.describe_group(aws_account_id, namespace, group_name)
        return group.add_member(member_name)

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
        self, aws_account_id: str, namespace: str, group_name: str, member_name: str
    ) -> QuicksightMembership:
        group = self.describe_group(aws_account_id, namespace, group_name)
        member = group.get_member(member_name)
        if member is None:
            raise ResourceNotFoundException(f"Member {member_name} not found")
        return member

    def describe_user(
        self, aws_account_id: str, namespace: str, user_name: str
    ) -> QuicksightUser:
        _id = _create_id(aws_account_id, namespace, user_name)
        if _id not in self.users:
            raise ResourceNotFoundException(f"User {user_name} not found")
        return self.users[_id]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_groups(self, aws_account_id: str, namespace: str) -> List[QuicksightGroup]:
        id_for_ns = _create_id(aws_account_id, namespace, _id="")
        return [
            group for _id, group in self.groups.items() if _id.startswith(id_for_ns)
        ]

    @paginate(pagination_model=PAGINATION_MODEL)
    def search_groups(
        self, aws_account_id: str, namespace: str, filters: List[Dict[str, str]]
    ) -> List[QuicksightGroup]:
        id_for_ns = _create_id(aws_account_id, namespace, _id="")
        filter_list = QuicksightSearchFilterFactory.validate_and_create_filter(
            model_type=QuicksightGroup, input=filters
        )
        return [
            group
            for _id, group in self.groups.items()
            if _id.startswith(id_for_ns) and filter_list.match(group)
        ]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_group_memberships(
        self, aws_account_id: str, namespace: str, group_name: str
    ) -> List[QuicksightMembership]:
        group = self.describe_group(aws_account_id, namespace, group_name)
        return group.list_members()

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_users(self, aws_account_id: str, namespace: str) -> List[QuicksightUser]:
        id_for_ns = _create_id(aws_account_id, namespace, _id="")
        return [user for _id, user in self.users.items() if _id.startswith(id_for_ns)]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_user_groups(
        self, aws_account_id: str, namespace: str, user_name: str
    ) -> List[QuicksightGroup]:
        id_for_ns = _create_id(aws_account_id, namespace, _id="")
        group_list: Dict[str, QuicksightGroup] = {}
        # Loop through all groups and check if the user is member.
        for id, group in self.groups.items():
            if group.get_member(user_name):
                group_list[id] = group
        return [group for _id, group in group_list.items() if _id.startswith(id_for_ns)]

    def register_user(
        self,
        identity_type: str,
        email: str,
        user_role: str,
        aws_account_id: str,
        namespace: str,
        user_name: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> QuicksightUser:
        """
        The following parameters are not yet implemented:
        IamArn, SessionName, CustomsPermissionsName, ExternalLoginFederationProviderType, CustomFederationProviderUrl, ExternalLoginId
        """
        user = QuicksightUser(
            account_id=aws_account_id,
            region=self.region_name,
            email=email,
            identity_type=identity_type,
            user_role=user_role,
            username=user_name,
        )
        _id = _create_id(aws_account_id, namespace, user_name)

        if tags:
            self.tagger.tag_resource(
                arn=user.arn,
                tags=tags,
            )

        self.users[_id] = user
        return user

    def update_user(
        self,
        aws_account_id: str,
        namespace: str,
        user_name: str,
        email: str,
        user_role: str,
    ) -> QuicksightUser:
        user = self.describe_user(aws_account_id, namespace, user_name)
        user.email = email
        user.user_role = user_role
        return user

    def update_group(
        self, aws_account_id: str, namespace: str, group_name: str, description: str
    ) -> QuicksightGroup:
        group = self.describe_group(aws_account_id, namespace, group_name)
        group.description = description
        return group

    def create_dashboard(
        self,
        aws_account_id: str,
        dashboard_id: str,
        name: str,
        parameters: Dict[str, Any],
        permissions: List[Dict[str, Any]],
        source_entity: Dict[str, Any],
        tags: List[Dict[str, str]],
        version_description: str,
        dashboard_publish_options: Dict[str, Any],
        theme_arn: str,
        definition: Dict[str, Any],
        validation_strategy: Dict[str, str],
        folder_arns: List[str],
        link_sharing_configuration: Dict[str, Any],
        link_entities: List[str],
    ) -> QuicksightDashboard:
        dashboard = QuicksightDashboard(
            account_id=aws_account_id,
            region=self.region_name,
            dashboard_id=dashboard_id,
            dashboard_publish_options=dashboard_publish_options,
            name=name,
            definition=definition,
            folder_arns=folder_arns,
            link_entities=link_entities,
            link_sharing_configuration=link_sharing_configuration,
            parameters=parameters,
            permissions=permissions,
            source_entity=source_entity,
            tags=tags,
            theme_arn=theme_arn,
            version_description=version_description,
            validation_strategy=validation_strategy,
        )

        if tags:
            self.tagger.tag_resource(
                arn=dashboard.arn,
                tags=tags,
            )

        self.dashboards[dashboard_id] = dashboard
        return dashboard

    def describe_dashboard(
        self,
        aws_account_id: str,
        dashboard_id: str,
        version_number: int,
        alias_name: str,
    ) -> QuicksightDashboard:
        dashboard = self.dashboards.get(dashboard_id)
        if not dashboard:
            raise ResourceNotFoundException(f"Dashboard {dashboard_id} not found")
        return dashboard

    def list_dashboards(self, aws_account_id: str) -> List[Dict[str, Any]]:
        dashboards = self.dashboards.values()
        dashboard_list: List[Dict[str, Any]] = []
        for dashboard in dashboards:
            d_dict = {
                "Arn": dashboard.arn,
                "DashboardId": dashboard.dashboard_id,
                "Name": dashboard.name,
                "CreatedTime": str(dashboard.created_time),
                "LastUpdatedTime": str(dashboard.last_updated_time),
                "PublishedVersionNumber": dashboard.version_number,
                "LastPublishedTime": str(dashboard.last_published_time),
            }
            dashboard_list.append(d_dict)

        return dashboard_list

    def describe_account_settings(
        self, aws_account_id: str
    ) -> QuicksightAccountSettings:
        return self.account_settings

    def update_account_settings(
        self,
        aws_account_id: str,
        default_namespace: str,
        notification_email: str,
        termination_protection_enabled: bool,
    ) -> None:
        if notification_email:
            self.account_settings.notification_email = notification_email
        if termination_protection_enabled:
            self.account_settings.termination_protection_enabled = (
                termination_protection_enabled
            )

    def update_public_sharing_settings(
        self, aws_account_id: str, public_sharing_enabled: bool
    ) -> None:
        self.account_settings.public_sharing_enabled = public_sharing_enabled

    def create_data_source(
        self,
        aws_account_id: str,
        data_source_id: str,
        name: str,
        data_source_type: str,
        data_source_parameters: Optional[Dict[str, Dict[str, Any]]] = None,
        vpc_connection_properties: Optional[Dict[str, Any]] = None,
        ssl_properties: Optional[Dict[str, Any]] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> QuickSightDataSource:
        data_source = QuickSightDataSource(
            account_id=aws_account_id,
            region=self.region_name,
            data_source_id=data_source_id,
            name=name,
            status="CREATION_SUCCESSFUL",
            data_source_type=data_source_type,
            data_source_parameters=data_source_parameters,
            ssl_properties=ssl_properties,
            vpc_connection_properties=vpc_connection_properties,
            tags=tags,
        )

        if tags:
            self.tagger.tag_resource(
                arn=data_source.arn,
                tags=tags,
            )

        self.data_sources[data_source_id] = data_source
        return data_source

    def update_data_source(
        self,
        aws_account_id: str,
        data_source_id: str,
        name: str,
        data_source_parameters: Optional[Dict[str, Any]] = None,
        vpc_connection_properties: Optional[Dict[str, Any]] = None,
        ssl_properties: Optional[Dict[str, Any]] = None,
    ) -> QuickSightDataSource:
        data_source = self.data_sources.get(data_source_id)

        if not data_source:
            raise ResourceNotFoundException(f"DataSource {data_source_id} Not Found")

        data_source.name = name
        data_source.data_source_parameters = data_source_parameters
        data_source.last_updated_time = datetime.now()
        data_source.status = "UPDATE_SUCCESSFUL"
        return data_source

    def delete_data_source(
        self, aws_account_id: str, data_source_id: str
    ) -> QuickSightDataSource:
        data_source = self.data_sources.pop(data_source_id, None)

        if data_source is None:
            raise ResourceNotFoundException(f"DataSource {data_source_id} Not Found")

        return data_source

    def describe_data_source(
        self, aws_account_id: str, data_source_id: str
    ) -> QuickSightDataSource:
        data_source = self.data_sources.get(data_source_id)

        if not data_source:
            raise ResourceNotFoundException(f"DataSource {data_source_id} Not Found")

        return data_source

    def list_data_sources(self, aws_account_id: str) -> List[Dict[str, Any]]:
        data_sources = self.data_sources.values()
        data_source_list: List[Dict[str, Any]] = []

        for data_source in data_sources:
            data_source_list.append(data_source.to_json())

        return data_source_list

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, str]]) -> None:
        self.tagger.tag_resource(
            arn=resource_arn,
            tags=tags,
        )

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(
            resource_arn,
            tag_keys,
        )

    def list_tags_for_resource(self, arn: str) -> List[Dict[str, str]]:
        tags = self.tagger.list_tags_for_resource(arn)
        return tags.get("Tags", [])


quicksight_backends = BackendDict(QuickSightBackend, "quicksight")
