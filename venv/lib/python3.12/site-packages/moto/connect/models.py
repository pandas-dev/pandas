"""ConnectBackend class with methods for supported APIs."""

import uuid
from datetime import datetime, timezone
from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService

from .exceptions import InvalidParameterException, ResourceNotFoundException

PAGINATION_MODEL = {
    "list_analytics_data_associations": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 1000,
        "unique_attribute": "DataSetId",
        "output_token": "next_token",
    },
    "list_instances": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 10,
        "unique_attribute": "Id",
        "output_token": "next_token",
    },
}


class Instance(BaseModel):
    def __init__(
        self,
        identity_management_type: str,
        inbound_calls_enabled: bool,
        outbound_calls_enabled: bool,
        account_id: str,
        region: str,
        instance_alias: Optional[str] = None,
        tags: Optional[dict[str, str]] = None,
    ) -> None:
        self.id = str(uuid.uuid4())
        self.arn = f"arn:aws:connect:{region}:{account_id}:instance/{self.id}"
        self.identity_management_type = identity_management_type
        self.instance_alias = instance_alias
        self.created_time = datetime.now(timezone.utc)
        self.service_role = f"arn:aws:iam::{account_id}:role/aws-service-role/connect.amazonaws.com/AWSServiceRoleForAmazonConnect"
        self.instance_status = "ACTIVE"
        self.inbound_calls_enabled = inbound_calls_enabled
        self.outbound_calls_enabled = outbound_calls_enabled
        self.instance_access_url = f"https://{instance_alias or self.id}.my.connect.aws"
        self.tags = tags or {}

    def to_dict(self) -> dict[str, Any]:
        """Return full instance details for DescribeInstance."""
        result: dict[str, Any] = {
            "Id": self.id,
            "Arn": self.arn,
            "IdentityManagementType": self.identity_management_type,
            "CreatedTime": self.created_time.isoformat(),
            "ServiceRole": self.service_role,
            "InstanceStatus": self.instance_status,
            "InboundCallsEnabled": self.inbound_calls_enabled,
            "OutboundCallsEnabled": self.outbound_calls_enabled,
            "InstanceAccessUrl": self.instance_access_url,
        }
        if self.instance_alias:
            result["InstanceAlias"] = self.instance_alias
        if self.tags:
            result["Tags"] = self.tags
        return result

    def to_summary_dict(self) -> dict[str, Any]:
        """Return summary for ListInstances."""
        result: dict[str, Any] = {
            "Id": self.id,
            "Arn": self.arn,
            "IdentityManagementType": self.identity_management_type,
            "CreatedTime": self.created_time.isoformat(),
            "ServiceRole": self.service_role,
            "InstanceStatus": self.instance_status,
            "InboundCallsEnabled": self.inbound_calls_enabled,
            "OutboundCallsEnabled": self.outbound_calls_enabled,
            "InstanceAccessUrl": self.instance_access_url,
        }
        if self.instance_alias:
            result["InstanceAlias"] = self.instance_alias
        return result


class AnalyticsDataAssociation(BaseModel):
    """Represents an analytics data association for a Connect instance.

    Storage model: One association per (instance_id, data_set_id) pair.
    A data_set_id can only be associated once per instance.
    The target_account_id specifies where the data is shared to.
    """

    def __init__(
        self,
        instance_id: str,
        data_set_id: str,
        target_account_id: str,
        source_account_id: str,
        region: str,
    ) -> None:
        self.instance_id = instance_id
        self.data_set_id = data_set_id
        self.target_account_id = target_account_id
        self.resource_share_id = str(uuid.uuid4())
        self.resource_share_arn = f"arn:aws:ram:{region}:{source_account_id}:resource-share/{self.resource_share_id}"
        self.resource_share_status = "ACTIVE"

    def to_dict(self) -> dict[str, str]:
        return {
            "DataSetId": self.data_set_id,
            "TargetAccountId": self.target_account_id,
            "ResourceShareId": self.resource_share_id,
            "ResourceShareArn": self.resource_share_arn,
            "ResourceShareStatus": self.resource_share_status,
        }


class ConnectBackend(BaseBackend):
    """Backend for Amazon Connect API."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.instances: dict[str, Instance] = {}
        self.analytics_data_associations: dict[
            str, dict[str, AnalyticsDataAssociation]
        ] = {}
        self.tagger = TaggingService()

    def tag_resource(self, resource_arn: str, tags: dict[str, str]) -> None:
        self.tagger.tag_resource(
            resource_arn, TaggingService.convert_dict_to_tags_input(tags)
        )

    def untag_resource(self, resource_arn: str, tag_keys: list[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def list_tags_for_resource(self, resource_arn: str) -> dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_arn)

    def _get_instance_or_raise(self, instance_id: str) -> Instance:
        """Get instance by ID or raise ResourceNotFoundException."""
        if instance_id not in self.instances:
            raise ResourceNotFoundException(f"Instance with id {instance_id} not found")
        return self.instances[instance_id]

    def create_instance(
        self,
        identity_management_type: str,
        inbound_calls_enabled: bool,
        outbound_calls_enabled: bool,
        instance_alias: Optional[str] = None,
        tags: Optional[dict[str, str]] = None,
    ) -> dict[str, str]:
        """Create a new Connect instance."""
        if not identity_management_type:
            raise InvalidParameterException(
                "IdentityManagementType is a required parameter"
            )

        instance = Instance(
            identity_management_type=identity_management_type,
            inbound_calls_enabled=inbound_calls_enabled,
            outbound_calls_enabled=outbound_calls_enabled,
            instance_alias=instance_alias,
            account_id=self.account_id,
            region=self.region_name,
            tags=tags,
        )
        self.instances[instance.id] = instance

        if tags:
            self.tag_resource(instance.arn, tags)

        return {
            "Id": instance.id,
            "Arn": instance.arn,
        }

    def describe_instance(self, instance_id: str) -> dict[str, Any]:
        """Describe a Connect instance."""
        if not instance_id:
            raise InvalidParameterException("InstanceId is a required parameter")

        instance = self._get_instance_or_raise(instance_id)
        return instance.to_dict()

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore[misc]
    def list_instances(
        self,
        max_results: Optional[int] = None,
        next_token: Optional[str] = None,
    ) -> list[dict[str, str]]:
        sorted_instances = sorted(
            self.instances.values(),
            key=lambda i: (i.created_time, i.id),
        )
        return [instance.to_summary_dict() for instance in sorted_instances]

    def delete_instance(self, instance_id: str) -> None:
        if not instance_id:
            raise InvalidParameterException("InstanceId is a required parameter")

        self._get_instance_or_raise(instance_id)

        del self.instances[instance_id]
        # Also clean up any analytics associations for this instance
        if instance_id in self.analytics_data_associations:
            del self.analytics_data_associations[instance_id]

    def associate_analytics_data_set(
        self,
        instance_id: str,
        data_set_id: str,
        target_account_id: Optional[str] = None,
    ) -> dict[str, str]:
        """Associate an analytics data set with a Connect instance.

        Enforces one association per (instance_id, data_set_id) pair.
        """
        if not instance_id:
            raise InvalidParameterException("InstanceId is a required parameter")
        if not data_set_id:
            raise InvalidParameterException("DataSetId is a required parameter")

        self._get_instance_or_raise(instance_id)

        if not target_account_id:
            target_account_id = self.account_id

        if instance_id in self.analytics_data_associations:
            if data_set_id in self.analytics_data_associations[instance_id]:
                raise InvalidParameterException(
                    f"Analytics data association for data set {data_set_id} already exists for instance {instance_id}"
                )

        if instance_id not in self.analytics_data_associations:
            self.analytics_data_associations[instance_id] = {}

        association = AnalyticsDataAssociation(
            instance_id=instance_id,
            data_set_id=data_set_id,
            target_account_id=target_account_id,
            source_account_id=self.account_id,
            region=self.region_name,
        )
        self.analytics_data_associations[instance_id][data_set_id] = association

        return {
            "DataSetId": association.data_set_id,
            "TargetAccountId": association.target_account_id,
            "ResourceShareId": association.resource_share_id,
            "ResourceShareArn": association.resource_share_arn,
            "ResourceShareStatus": association.resource_share_status,
        }

    def disassociate_analytics_data_set(
        self,
        instance_id: str,
        data_set_id: str,
    ) -> None:
        """Disassociate an analytics data set from a Connect instance."""
        if not instance_id:
            raise InvalidParameterException("InstanceId is a required parameter")
        if not data_set_id:
            raise InvalidParameterException("DataSetId is a required parameter")

        self._get_instance_or_raise(instance_id)

        if instance_id not in self.analytics_data_associations:
            raise ResourceNotFoundException(
                f"Analytics data association for data set {data_set_id} not found for instance {instance_id}"
            )

        if data_set_id not in self.analytics_data_associations[instance_id]:
            raise ResourceNotFoundException(
                f"Analytics data association for data set {data_set_id} not found for instance {instance_id}"
            )

        del self.analytics_data_associations[instance_id][data_set_id]

        if not self.analytics_data_associations[instance_id]:
            del self.analytics_data_associations[instance_id]

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore[misc]
    def list_analytics_data_associations(
        self,
        instance_id: str,
        data_set_id: Optional[str] = None,
        max_results: Optional[int] = None,
        next_token: Optional[str] = None,
    ) -> list[dict[str, str]]:
        """List analytics data associations for a Connect instance."""
        if not instance_id:
            raise InvalidParameterException("InstanceId is a required parameter")

        self._get_instance_or_raise(instance_id)

        if instance_id not in self.analytics_data_associations:
            return []

        associations = list(self.analytics_data_associations[instance_id].values())

        if data_set_id:
            associations = [a for a in associations if a.data_set_id == data_set_id]

        associations.sort(key=lambda a: a.data_set_id)

        return [a.to_dict() for a in associations]


connect_backends = BackendDict(ConnectBackend, "connect")
