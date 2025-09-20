"""Handles incoming quicksight requests, invokes methods, returns responses."""

import json
from urllib.parse import unquote

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import QuickSightBackend, quicksight_backends


class QuickSightResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="quicksight")

    @property
    def quicksight_backend(self) -> QuickSightBackend:
        """Return backend instance specific for this region."""
        return quicksight_backends[self.current_account][self.region]

    def create_data_set(self) -> str:
        params = json.loads(self.body)
        data_set_id = params.get("DataSetId")
        name = params.get("Name")
        tags = self._get_param("Tags")

        data_set = self.quicksight_backend.create_data_set(data_set_id, name, tags)
        return json.dumps(data_set.to_json())

    def create_group(self) -> str:
        params = json.loads(self.body)
        group_name = params.get("GroupName")
        description = params.get("Description")
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        group = self.quicksight_backend.create_group(
            group_name=group_name,
            description=description,
            aws_account_id=aws_account_id,
            namespace=namespace,
        )
        return json.dumps(dict(Group=group.to_json()))

    def create_group_membership(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        group_name = unquote(self._get_param("GroupName"))
        member_name = unquote(self._get_param("MemberName"))
        member = self.quicksight_backend.create_group_membership(
            aws_account_id, namespace, group_name, member_name
        )
        return json.dumps({"GroupMember": member.to_json()})

    def create_ingestion(self) -> str:
        data_set_id = self._get_param("DataSetId")
        ingestion_id = self._get_param("IngestionId")
        ingestion = self.quicksight_backend.create_ingestion(data_set_id, ingestion_id)
        return json.dumps(ingestion.to_json())

    def describe_group_membership(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        group_name = unquote(self._get_param("GroupName"))
        member_name = unquote(self._get_param("MemberName"))
        member = self.quicksight_backend.describe_group_membership(
            aws_account_id, namespace, group_name, member_name
        )
        return json.dumps({"GroupMember": member.to_json()})

    def list_groups(self) -> str:
        max_results = self._get_int_param("max-results")
        next_token = self._get_param("next-token")
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        groups, next_token = self.quicksight_backend.list_groups(
            aws_account_id, namespace, max_results=max_results, next_token=next_token
        )
        return json.dumps(
            {"NextToken": next_token, "GroupList": [g.to_json() for g in groups]}
        )

    def list_group_memberships(self) -> str:
        max_results = self._get_int_param("max-results")
        next_token = self._get_param("next-token")
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        group_name = unquote(self._get_param("GroupName"))
        members, next_token = self.quicksight_backend.list_group_memberships(
            aws_account_id,
            namespace,
            group_name,
            max_results=max_results,
            next_token=next_token,
        )
        return json.dumps(
            {"NextToken": next_token, "GroupMemberList": [m.to_json() for m in members]}
        )

    def list_users(self) -> str:
        max_results = self._get_int_param("max-results")
        next_token = self._get_param("next-token")
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        users, next_token = self.quicksight_backend.list_users(
            aws_account_id, namespace, max_results=max_results, next_token=next_token
        )
        return json.dumps(
            {"NextToken": next_token, "UserList": [u.to_json() for u in users]}
        )

    def list_user_groups(self) -> str:
        max_results = self._get_int_param("max-results")
        next_token = self._get_param("next-token")
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        user_name = unquote(self._get_param("UserName"))
        groups, next_token = self.quicksight_backend.list_user_groups(
            aws_account_id,
            namespace,
            user_name,
            max_results=max_results,
            next_token=next_token,
        )
        return json.dumps(
            {"NextToken": next_token, "GroupList": [g.to_json() for g in groups]}
        )

    def register_user(self) -> str:
        params = json.loads(self.body)
        identity_type = params.get("IdentityType")
        email = params.get("Email")
        user_role = params.get("UserRole")
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        user_name = params.get("UserName")
        tags = params.get("Tags")
        user = self.quicksight_backend.register_user(
            identity_type=identity_type,
            email=email,
            user_role=user_role,
            aws_account_id=aws_account_id,
            namespace=namespace,
            user_name=user_name,
            tags=tags,
        )
        return json.dumps(dict(User=user.to_json(), UserInvitationUrl="TBD"))

    def update_user(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        user_name = unquote(self._get_param("UserName"))
        body = json.loads(self.body)
        email = body.get("Email", None)
        user_role = body.get("Role", None)

        user = self.quicksight_backend.update_user(
            aws_account_id, namespace, user_name, email, user_role
        )
        return json.dumps(dict(User=user.to_json()))

    def describe_group(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        group_name = unquote(self._get_param("GroupName"))

        group = self.quicksight_backend.describe_group(
            aws_account_id, namespace, group_name
        )
        return json.dumps(dict(Group=group.to_json()))

    def describe_user(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        user_name = unquote(self._get_param("UserName"))

        user = self.quicksight_backend.describe_user(
            aws_account_id, namespace, user_name
        )
        return json.dumps(dict(User=user.to_json()))

    def delete_group(self) -> TYPE_RESPONSE:
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        group_name = unquote(self._get_param("GroupName"))

        self.quicksight_backend.delete_group(aws_account_id, namespace, group_name)
        return 204, {"status": 204}, json.dumps({"Status": 204})

    def delete_user(self) -> TYPE_RESPONSE:
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        user_name = unquote(self._get_param("UserName"))

        self.quicksight_backend.delete_user(aws_account_id, namespace, user_name)
        return 204, {"status": 204}, json.dumps({"Status": 204})

    def update_group(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        group_name = unquote(self._get_param("GroupName"))
        description = json.loads(self.body).get("Description")

        group = self.quicksight_backend.update_group(
            aws_account_id, namespace, group_name, description
        )
        return json.dumps(dict(Group=group.to_json()))

    def search_groups(self) -> str:
        max_results = self._get_int_param("max-results")
        next_token = self._get_param("next-token")
        aws_account_id = self._get_param("AwsAccountId")
        namespace = self._get_param("Namespace")
        body = json.loads(self.body)

        groups, next_token = self.quicksight_backend.search_groups(
            aws_account_id,
            namespace,
            body.get("Filters", None),
            max_results=max_results,
            next_token=next_token,
        )
        return json.dumps(
            {"NextToken": next_token, "GroupList": [g.to_json() for g in groups]}
        )

    def create_dashboard(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        dashboard_id = self._get_param("DashboardId")
        dashboard_publish_options = self._get_param("DashboardPublishOptions")
        definition = self._get_param("Definition")
        folder_arns = self._get_param("FolderArns")
        link_entities = self._get_param("LinkEntities")
        link_sharing_configuration = self._get_param("LinkSharingConfiguration")
        name = self._get_param("Name")
        parameters = self._get_param("Parameters")
        permissions = self._get_param("Permissions")
        source_entity = self._get_param("SourceEntity")
        tags = self._get_param("Tags")
        theme_arn = self._get_param("ThemeArn")
        validation_strategy = self._get_param("ValidationStrategy")
        version_description = self._get_param("VersionDescription")

        dashboard = self.quicksight_backend.create_dashboard(
            aws_account_id=aws_account_id,
            dashboard_id=dashboard_id,
            name=name,
            parameters=parameters,
            permissions=permissions,
            source_entity=source_entity,
            tags=tags,
            version_description=version_description,
            dashboard_publish_options=dashboard_publish_options,
            theme_arn=theme_arn,
            definition=definition,
            validation_strategy=validation_strategy,
            folder_arns=folder_arns,
            link_sharing_configuration=link_sharing_configuration,
            link_entities=link_entities,
        )
        return json.dumps(
            dict(
                Arn=dashboard.arn,
                VersionArn=dashboard.version_number,
                DashboardId=dashboard.dashboard_id,
                CreationStatus=dashboard.status,
            )
        )

    def describe_dashboard(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        dashboard_id = self._get_param("DashboardId")
        version_number = self._get_param("VersionNumber")
        alias_name = self._get_param("AliasName")
        dashboard = self.quicksight_backend.describe_dashboard(
            aws_account_id=aws_account_id,
            dashboard_id=dashboard_id,
            version_number=version_number,
            alias_name=alias_name,
        )
        return json.dumps(
            dict(Dashboard=dashboard.to_dict(), Status=200, RequestId="request_id")
        )

    def list_dashboards(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        next_token = self._get_param("NextToken")
        dashboard_summary_list = self.quicksight_backend.list_dashboards(
            aws_account_id=aws_account_id,
        )
        return json.dumps(
            dict(
                DashboardSummaryList=dashboard_summary_list,
                Next_token=next_token,
                Status=200,
            )
        )

    def describe_account_settings(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        settings = self.quicksight_backend.describe_account_settings(
            aws_account_id=aws_account_id,
        )
        resp = {
            "AccountName": settings.account_name,
            "Edition": settings.edition,
            "DefaultNamespace": settings.default_namespace,
            "NotificationEmail": settings.notification_email,
            "PublicSharingEnabled": settings.public_sharing_enabled,
            "TerminationProtectionEnabled": settings.termination_protection_enabled,
        }

        return json.dumps(dict(AccountSettings=resp, Status=200))

    def update_account_settings(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        default_namespace = self._get_param("DefaultNamespace")
        notification_email = self._get_param("NotificationEmail")
        termination_protection_enabled = self._get_param("TerminationProtectionEnabled")
        self.quicksight_backend.update_account_settings(
            aws_account_id=aws_account_id,
            default_namespace=default_namespace,
            notification_email=notification_email,
            termination_protection_enabled=termination_protection_enabled,
        )
        return json.dumps(dict(Status=200))

    def update_public_sharing_settings(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        public_sharing_enabled = self._get_param("PublicSharingEnabled")
        self.quicksight_backend.update_public_sharing_settings(
            aws_account_id=aws_account_id,
            public_sharing_enabled=public_sharing_enabled,
        )
        return json.dumps(dict(Status=200))

    def create_data_source(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        data_source_id = self._get_param("DataSourceId")
        name = self._get_param("Name")
        data_source_type = self._get_param("Type")
        data_source_parameters = self._get_param("DataSourceParameters")
        ssl_properties = self._get_param("SslProperties")
        vpc_connection_properties = self._get_param("VpcConnectionProperties")
        tags = self._get_param("Tags")

        data_source = self.quicksight_backend.create_data_source(
            aws_account_id=aws_account_id,
            data_source_id=data_source_id,
            name=name,
            data_source_type=data_source_type,
            data_source_parameters=data_source_parameters,
            ssl_properties=ssl_properties,
            vpc_connection_properties=vpc_connection_properties,
            tags=tags,
        )

        return json.dumps(
            {
                "Arn": data_source.arn,
                "DataSourceId": data_source.data_source_id,
                "CreationStatus": data_source.status,
                "RequestId": "request_id",
                "Status": 200,
            }
        )

    def delete_data_source(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        data_source_id = self._get_param("DataSourceId")

        resp = self.quicksight_backend.delete_data_source(
            aws_account_id=aws_account_id, data_source_id=data_source_id
        )

        return json.dumps(
            {
                "Arn": resp.arn,
                "DataSourceId": resp.data_source_id,
                "RequestId": "request_id",
                "Status": 200,
            }
        )

    def describe_data_source(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        data_source_id = self._get_param("DataSourceId")

        data_source = self.quicksight_backend.describe_data_source(
            aws_account_id=aws_account_id, data_source_id=data_source_id
        )

        return json.dumps(
            {
                "DataSource": data_source.to_json(),
                "RequestId": "request_id",
                "Status": 200,
            }
        )

    def list_data_sources(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        next_token = self._get_param("NextToken")

        data_sources = self.quicksight_backend.list_data_sources(
            aws_account_id=aws_account_id,
        )

        return json.dumps(
            {
                "DataSources": data_sources,
                "NextToken": next_token,
                "RequestId": "request_id",
                "Status": 200,
            }
        )

    def update_data_source(self) -> str:
        aws_account_id = self._get_param("AwsAccountId")
        data_source_id = self._get_param("DataSourceId")
        name = self._get_param("Name")
        data_source_parameters = self._get_param("DataSourceParameters")

        data_source = self.quicksight_backend.update_data_source(
            aws_account_id=aws_account_id,
            data_source_id=data_source_id,
            name=name,
            data_source_parameters=data_source_parameters,
        )

        return json.dumps(
            {
                "Arn": data_source.arn,
                "DataSourceId": data_source.data_source_id,
                "UpdateStatus": data_source.status,
                "RequestId": "request_id",
                "Status": 200,
            }
        )

    def tag_resource(self) -> str:
        resource_arn = unquote(self._get_param("ResourceArn"))
        tags = self._get_param("Tags")

        self.quicksight_backend.tag_resource(resource_arn, tags)

        return json.dumps(dict(RequestId="request_id", Status=200))

    def untag_resource(self) -> str:
        resource_arn = unquote(self._get_param("ResourceArn"))
        tag_keys = self.__dict__["data"]["keys"]

        self.quicksight_backend.untag_resource(resource_arn, tag_keys)

        return json.dumps(dict(RequestId="request_id", Status=200))

    def list_tags_for_resource(self) -> str:
        resource_arn = unquote(self._get_param("ResourceArn"))
        tags = self.quicksight_backend.list_tags_for_resource(arn=resource_arn)

        return json.dumps(dict(Tags=tags, RequestId="request_id", Status=200))
