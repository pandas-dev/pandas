"""PrometheusServiceBackend class with methods for supported APIs."""

from typing import Any, Callable, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.moto_api._internal import mock_random
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import RuleGroupNamespaceNotFound, WorkspaceNotFound
from .utils import PAGINATION_MODEL


class RuleGroupNamespace(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        workspace_id: str,
        name: str,
        data: str,
        tag_fn: Callable[[str], Dict[str, str]],
    ):
        self.name = name
        self.data = data
        self.tag_fn = tag_fn
        self.arn = f"arn:{get_partition(region)}:aps:{region}:{account_id}:rulegroupsnamespace/{workspace_id}/{self.name}"
        self.created_at = unix_time()
        self.modified_at = self.created_at

    def update(self, new_data: str) -> None:
        self.data = new_data
        self.modified_at = unix_time()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "arn": self.arn,
            "status": {"statusCode": "ACTIVE"},
            "createdAt": self.created_at,
            "modifiedAt": self.modified_at,
            "data": self.data,
            "tags": self.tag_fn(self.arn),
        }


class Workspace(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        alias: str,
        tag_fn: Callable[[str], Dict[str, str]],
    ):
        self.alias = alias
        self.workspace_id = f"ws-{mock_random.uuid4()}"
        self.arn = f"arn:{get_partition(region)}:aps:{region}:{account_id}:workspace/{self.workspace_id}"
        self.endpoint = f"https://aps-workspaces.{region}.amazonaws.com/workspaces/{self.workspace_id}/"
        self.status = {"statusCode": "ACTIVE"}
        self.created_at = unix_time()
        self.tag_fn = tag_fn
        self.rule_group_namespaces: Dict[str, RuleGroupNamespace] = dict()
        self.logging_config: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "alias": self.alias,
            "arn": self.arn,
            "workspaceId": self.workspace_id,
            "status": self.status,
            "createdAt": self.created_at,
            "prometheusEndpoint": self.endpoint,
            "tags": self.tag_fn(self.arn),
        }


class PrometheusServiceBackend(BaseBackend):
    """Implementation of PrometheusService APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.workspaces: Dict[str, Workspace] = dict()
        self.tagger = TaggingService()

    def create_workspace(self, alias: str, tags: Dict[str, str]) -> Workspace:
        """
        The ClientToken-parameter is not yet implemented
        """
        workspace = Workspace(
            self.account_id,
            self.region_name,
            alias=alias,
            tag_fn=self.list_tags_for_resource,
        )
        self.workspaces[workspace.workspace_id] = workspace
        self.tag_resource(workspace.arn, tags)
        return workspace

    def describe_workspace(self, workspace_id: str) -> Workspace:
        if workspace_id not in self.workspaces:
            raise WorkspaceNotFound(workspace_id)
        return self.workspaces[workspace_id]

    def list_tags_for_resource(self, resource_arn: str) -> Dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_arn)

    def update_workspace_alias(self, alias: str, workspace_id: str) -> None:
        """
        The ClientToken-parameter is not yet implemented
        """
        self.workspaces[workspace_id].alias = alias

    def delete_workspace(self, workspace_id: str) -> None:
        """
        The ClientToken-parameter is not yet implemented
        """
        self.workspaces.pop(workspace_id, None)

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_workspaces(self, alias: str) -> List[Workspace]:
        if alias:
            return [w for w in self.workspaces.values() if w.alias == alias]
        return list(self.workspaces.values())

    def tag_resource(self, resource_arn: str, tags: Dict[str, str]) -> None:
        tag_list = self.tagger.convert_dict_to_tags_input(tags)
        self.tagger.tag_resource(resource_arn, tag_list)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def create_rule_groups_namespace(
        self, data: str, name: str, tags: Dict[str, str], workspace_id: str
    ) -> RuleGroupNamespace:
        """
        The ClientToken-parameter is not yet implemented
        """
        workspace = self.describe_workspace(workspace_id)
        group = RuleGroupNamespace(
            account_id=self.account_id,
            region=self.region_name,
            workspace_id=workspace_id,
            name=name,
            data=data,
            tag_fn=self.list_tags_for_resource,
        )
        workspace.rule_group_namespaces[name] = group
        self.tag_resource(group.arn, tags)
        return group

    def delete_rule_groups_namespace(self, name: str, workspace_id: str) -> None:
        """
        The ClientToken-parameter is not yet implemented
        """
        ws = self.describe_workspace(workspace_id)
        ws.rule_group_namespaces.pop(name, None)

    def describe_rule_groups_namespace(
        self, name: str, workspace_id: str
    ) -> RuleGroupNamespace:
        ws = self.describe_workspace(workspace_id)
        if name not in ws.rule_group_namespaces:
            raise RuleGroupNamespaceNotFound(name=name)
        return ws.rule_group_namespaces[name]

    def put_rule_groups_namespace(
        self, data: str, name: str, workspace_id: str
    ) -> RuleGroupNamespace:
        """
        The ClientToken-parameter is not yet implemented
        """
        ns = self.describe_rule_groups_namespace(name=name, workspace_id=workspace_id)
        ns.update(data)
        return ns

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_rule_groups_namespaces(
        self, name: str, workspace_id: str
    ) -> List[RuleGroupNamespace]:
        ws = self.describe_workspace(workspace_id)
        if name:
            return [
                ns
                for ns_name, ns in ws.rule_group_namespaces.items()
                if ns_name.startswith(name)
            ]
        return list(ws.rule_group_namespaces.values())

    def create_logging_configuration(
        self, workspace_id: str, log_group_arn: str
    ) -> Dict[str, str]:
        ws = self.describe_workspace(workspace_id)
        ws.logging_config = {
            "logGroupArn": log_group_arn,
            "createdAt": unix_time(),
            "status": {"statusCode": "ACTIVE"},
            "workspace": workspace_id,
        }
        return ws.logging_config["status"]

    def describe_logging_configuration(self, workspace_id: str) -> Dict[str, Any]:
        ws = self.describe_workspace(workspace_id)
        if ws.logging_config is None:
            return {}
        return ws.logging_config

    def delete_logging_configuration(self, workspace_id: str) -> None:
        ws = self.describe_workspace(workspace_id)
        ws.logging_config = None

    def update_logging_configuration(
        self, workspace_id: str, log_group_arn: str
    ) -> Dict[str, str]:
        ws = self.describe_workspace(workspace_id)
        ws.logging_config["logGroupArn"] = log_group_arn  # type: ignore[index]
        ws.logging_config["modifiedAt"] = unix_time()  # type: ignore[index]
        return ws.logging_config["status"]  # type: ignore[index]


amp_backends = BackendDict(PrometheusServiceBackend, "amp")
