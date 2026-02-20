"""Handles incoming connect requests, invokes methods, returns responses."""

import json
from typing import Any, Optional
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import ConnectBackend, connect_backends


class ConnectResponse(BaseResponse):
    """Handler for Connect requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="connect")

    @property
    def connect_backend(self) -> ConnectBackend:
        """Return backend instance specific for this region."""
        return connect_backends[self.current_account][self.region]

    def _get_instance_id(self) -> str:
        """Extract instance_id from request path params."""
        instance_id = self._get_param("InstanceId")
        return unquote(instance_id) if instance_id else ""

    def associate_analytics_data_set(self) -> str:
        instance_id = self._get_instance_id()
        params = json.loads(self.body) if self.body else {}
        data_set_id = str(params["DataSetId"])
        target_account_id = params.get("TargetAccountId")

        result = self.connect_backend.associate_analytics_data_set(
            instance_id=instance_id,
            data_set_id=data_set_id,
            target_account_id=target_account_id,
        )

        return json.dumps(result)

    def disassociate_analytics_data_set(self) -> str:
        instance_id = self._get_instance_id()
        params = json.loads(self.body) if self.body else {}
        data_set_id = str(params["DataSetId"])

        self.connect_backend.disassociate_analytics_data_set(
            instance_id=instance_id,
            data_set_id=data_set_id,
        )

        return "{}"

    def list_analytics_data_associations(self) -> str:
        instance_id = self._get_instance_id()
        data_set_id = self._get_param("DataSetId")
        max_results = self._get_int_param("maxResults")
        next_token = self._get_param("nextToken")

        results: list[dict[str, str]]
        token: Optional[str]
        results, token = self.connect_backend.list_analytics_data_associations(
            instance_id=instance_id,
            data_set_id=data_set_id,
            max_results=max_results,
            next_token=next_token,
        )

        response: dict[str, Any] = {"Results": results}
        if token:
            response["NextToken"] = token

        return json.dumps(response)

    def create_instance(self) -> str:
        params = json.loads(self.body) if self.body else {}
        identity_management_type = str(params["IdentityManagementType"])
        instance_alias = params.get("InstanceAlias")
        inbound_calls_enabled = params.get("InboundCallsEnabled", False)
        outbound_calls_enabled = params.get("OutboundCallsEnabled", False)
        tags = params.get("Tags")

        result = self.connect_backend.create_instance(
            identity_management_type=identity_management_type,
            instance_alias=instance_alias,
            inbound_calls_enabled=inbound_calls_enabled,
            outbound_calls_enabled=outbound_calls_enabled,
            tags=tags,
        )

        return json.dumps(result)

    def describe_instance(self) -> str:
        instance_id = self._get_instance_id()

        instance = self.connect_backend.describe_instance(instance_id=instance_id)

        return json.dumps({"Instance": instance})

    def list_instances(self) -> str:
        max_results = self._get_int_param("maxResults")
        next_token = self._get_param("nextToken")

        results: list[dict[str, Any]]
        token: Optional[str]
        results, token = self.connect_backend.list_instances(
            max_results=max_results,
            next_token=next_token,
        )

        response: dict[str, Any] = {"InstanceSummaryList": results}
        if token:
            response["NextToken"] = token

        return json.dumps(response)

    def delete_instance(self) -> str:
        instance_id = self._get_instance_id()

        self.connect_backend.delete_instance(instance_id=instance_id)

        return "{}"

    def _get_resource_arn(self) -> str:
        """Extract resourceArn from request path params."""
        resource_arn = self._get_param("resourceArn")
        return unquote(resource_arn) if resource_arn else ""

    def tag_resource(self) -> str:
        resource_arn = self._get_resource_arn()
        params = json.loads(self.body) if self.body else {}
        tags = params.get("tags", {})

        self.connect_backend.tag_resource(resource_arn=resource_arn, tags=tags)

        return "{}"

    def untag_resource(self) -> str:
        resource_arn = self._get_resource_arn()
        tag_keys = self.querystring.get("tagKeys", [])

        self.connect_backend.untag_resource(
            resource_arn=resource_arn, tag_keys=tag_keys
        )

        return "{}"

    def list_tags_for_resource(self) -> str:
        resource_arn = self._get_resource_arn()

        tags = self.connect_backend.list_tags_for_resource(resource_arn=resource_arn)

        return json.dumps({"tags": tags})
