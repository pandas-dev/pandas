"""Handles incoming amp requests, invokes methods, returns responses."""

import json
from typing import Any
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import PrometheusServiceBackend, amp_backends


class PrometheusServiceResponse(BaseResponse):
    """Handler for PrometheusService requests and responses."""

    def tags(self, request: Any, full_url: str, headers: Any) -> str:  # type: ignore[return]
        self.setup_class(request, full_url, headers)
        if request.method == "GET":
            return self.list_tags_for_resource()
        if request.method == "POST":
            return self.tag_resource()
        if request.method == "DELETE":
            return self.untag_resource()

    def __init__(self) -> None:
        super().__init__(service_name="amp")

    @property
    def amp_backend(self) -> PrometheusServiceBackend:
        """Return backend instance specific for this region."""
        return amp_backends[self.current_account][self.region]

    def create_workspace(self) -> str:
        params = json.loads(self.body)
        alias = params.get("alias")
        tags = params.get("tags")
        workspace = self.amp_backend.create_workspace(alias=alias, tags=tags)
        return json.dumps(dict(workspace.to_dict()))

    def describe_workspace(self) -> str:
        workspace_id = self.path.split("/")[-1]
        workspace = self.amp_backend.describe_workspace(workspace_id=workspace_id)
        return json.dumps(dict(workspace=workspace.to_dict()))

    def list_tags_for_resource(self) -> str:
        resource_arn = unquote(self.path).split("tags/")[-1]
        tags = self.amp_backend.list_tags_for_resource(resource_arn=resource_arn)
        return json.dumps(dict(tags=tags))

    def update_workspace_alias(self) -> str:
        params = json.loads(self.body)
        alias = params.get("alias")
        workspace_id = self.path.split("/")[-2]
        self.amp_backend.update_workspace_alias(alias=alias, workspace_id=workspace_id)
        return json.dumps(dict())

    def delete_workspace(self) -> str:
        workspace_id = self.path.split("/")[-1]
        self.amp_backend.delete_workspace(workspace_id=workspace_id)
        return json.dumps(dict())

    def list_workspaces(self) -> str:
        alias = self._get_param("alias")
        max_results = self._get_int_param("maxResults")
        next_token = self._get_param("nextToken")
        workspaces, next_token = self.amp_backend.list_workspaces(
            alias, max_results=max_results, next_token=next_token
        )
        return json.dumps(
            {"nextToken": next_token, "workspaces": [w.to_dict() for w in workspaces]}
        )

    def tag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = unquote(self.path).split("tags/")[-1]
        tags = params.get("tags")
        self.amp_backend.tag_resource(resource_arn=resource_arn, tags=tags)
        return json.dumps(dict())

    def untag_resource(self) -> str:
        resource_arn = unquote(self.path).split("tags/")[-1]
        tag_keys = self.querystring.get("tagKeys", [])
        self.amp_backend.untag_resource(resource_arn=resource_arn, tag_keys=tag_keys)
        return json.dumps(dict())

    def create_rule_groups_namespace(self) -> str:
        params = json.loads(self.body)
        data = params.get("data")
        name = params.get("name")
        tags = params.get("tags")
        workspace_id = unquote(self.path).split("/")[-2]
        rule_group_namespace = self.amp_backend.create_rule_groups_namespace(
            data=data,
            name=name,
            tags=tags,
            workspace_id=workspace_id,
        )
        return json.dumps(rule_group_namespace.to_dict())

    def delete_rule_groups_namespace(self) -> str:
        name = unquote(self.path).split("/")[-1]
        workspace_id = unquote(self.path).split("/")[-3]
        self.amp_backend.delete_rule_groups_namespace(
            name=name,
            workspace_id=workspace_id,
        )
        return json.dumps(dict())

    def describe_rule_groups_namespace(self) -> str:
        name = unquote(self.path).split("/")[-1]
        workspace_id = unquote(self.path).split("/")[-3]
        ns = self.amp_backend.describe_rule_groups_namespace(
            name=name, workspace_id=workspace_id
        )
        return json.dumps(dict(ruleGroupsNamespace=ns.to_dict()))

    def put_rule_groups_namespace(self) -> str:
        params = json.loads(self.body)
        data = params.get("data")
        name = unquote(self.path).split("/")[-1]
        workspace_id = unquote(self.path).split("/")[-3]
        ns = self.amp_backend.put_rule_groups_namespace(
            data=data,
            name=name,
            workspace_id=workspace_id,
        )
        return json.dumps(ns.to_dict())

    def list_rule_groups_namespaces(self) -> str:
        max_results = self._get_int_param("maxResults")
        next_token = self._get_param("nextToken")
        name = self._get_param("name")
        workspace_id = unquote(self.path).split("/")[-2]
        namespaces, next_token = self.amp_backend.list_rule_groups_namespaces(
            max_results=max_results,
            name=name,
            next_token=next_token,
            workspace_id=workspace_id,
        )
        return json.dumps(
            dict(
                nextToken=next_token,
                ruleGroupsNamespaces=[ns.to_dict() for ns in namespaces],
            )
        )

    def create_logging_configuration(self) -> str:
        workspace_id = unquote(self.path).split("/")[-2]
        log_group_arn = self._get_param("logGroupArn")
        status = self.amp_backend.create_logging_configuration(
            workspace_id=workspace_id,
            log_group_arn=log_group_arn,
        )
        return json.dumps({"status": status})

    def describe_logging_configuration(self) -> str:
        workspace_id = unquote(self.path).split("/")[-2]
        config = self.amp_backend.describe_logging_configuration(
            workspace_id=workspace_id
        )
        return json.dumps({"loggingConfiguration": config})

    def update_logging_configuration(self) -> str:
        workspace_id = unquote(self.path).split("/")[-2]
        log_group_arn = self._get_param("logGroupArn")
        status = self.amp_backend.update_logging_configuration(
            workspace_id=workspace_id, log_group_arn=log_group_arn
        )
        return json.dumps({"status": status})

    def delete_logging_configuration(self) -> str:
        workspace_id = unquote(self.path).split("/")[-2]
        self.amp_backend.delete_logging_configuration(workspace_id=workspace_id)
        return "{}"
