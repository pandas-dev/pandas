import json
from typing import Any, Dict, List
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import Inspector2Backend, inspector2_backends


class Inspector2Response(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="inspector2")

    @property
    def inspector2_backend(self) -> Inspector2Backend:
        return inspector2_backends[self.current_account][self.region]

    def create_filter(self) -> str:
        action = self._get_param("action")
        description = self._get_param("description")
        filter_criteria = self._get_param("filterCriteria")
        name = self._get_param("name")
        reason = self._get_param("reason")
        tags = self._get_param("tags")
        arn = self.inspector2_backend.create_filter(
            action=action,
            description=description,
            filter_criteria=filter_criteria,
            name=name,
            reason=reason,
            tags=tags,
        )
        return json.dumps(dict(arn=arn))

    def delete_filter(self) -> str:
        arn = self._get_param("arn")
        self.inspector2_backend.delete_filter(arn=arn)
        return json.dumps(dict(arn=arn))

    def list_filters(self) -> str:
        action = self._get_param("action")
        arns = self._get_param("arns")
        filters = self.inspector2_backend.list_filters(action=action, arns=arns)
        return json.dumps({"filters": [f.to_json() for f in filters]})

    def list_findings(self) -> str:
        filter_criteria = self._get_param("filterCriteria")
        max_results = self._get_param("maxResults")
        next_token = self._get_param("nextToken")
        sort_criteria = self._get_param("sortCriteria")
        findings = self.inspector2_backend.list_findings(
            filter_criteria=filter_criteria,
            max_results=max_results,
            next_token=next_token,
            sort_criteria=sort_criteria,
        )
        return json.dumps(dict(findings=findings))

    def list_delegated_admin_accounts(self) -> str:
        accounts = self.inspector2_backend.list_delegated_admin_accounts()
        return json.dumps(
            {
                "delegatedAdminAccounts": [
                    {"accountId": key, "status": val} for key, val in accounts.items()
                ]
            }
        )

    def enable_delegated_admin_account(self) -> str:
        account_id = self._get_param("delegatedAdminAccountId")
        self.inspector2_backend.enable_delegated_admin_account(account_id)
        return json.dumps({"delegatedAdminAccountId": account_id})

    def disable_delegated_admin_account(self) -> str:
        account_id = self._get_param("delegatedAdminAccountId")
        self.inspector2_backend.disable_delegated_admin_account(account_id)
        return json.dumps({"delegatedAdminAccountId": account_id})

    def describe_organization_configuration(self) -> str:
        config = self.inspector2_backend.describe_organization_configuration()
        return json.dumps(config)

    def update_organization_configuration(self) -> str:
        auto_enable = self._get_param("autoEnable")
        config = self.inspector2_backend.update_organization_configuration(auto_enable)
        return json.dumps(config)

    def enable(self) -> str:
        account_ids = self._get_param("accountIds")
        resource_types = self._get_param("resourceTypes")
        accounts = self.inspector2_backend.enable(account_ids, resource_types)
        failed: List[Dict[str, Any]] = []
        return json.dumps({"accounts": accounts, "failedAccounts": failed})

    def disable(self) -> str:
        account_ids = self._get_param("accountIds")
        resource_types = self._get_param("resourceTypes")
        accounts = self.inspector2_backend.disable(account_ids, resource_types)
        failed: List[Dict[str, Any]] = []
        return json.dumps({"accounts": accounts, "failedAccounts": failed})

    def batch_get_account_status(self) -> str:
        account_ids = self._get_param("accountIds")
        accounts = self.inspector2_backend.batch_get_account_status(account_ids)
        failed: List[Dict[str, Any]] = []
        return json.dumps({"accounts": accounts, "failedAccounts": failed})

    def list_members(self) -> str:
        members = self.inspector2_backend.list_members()
        return json.dumps({"members": [m.to_json() for m in members]})

    def associate_member(self) -> str:
        account_id = self._get_param("accountId")
        self.inspector2_backend.associate_member(account_id)
        return json.dumps({"accountId": account_id})

    def disassociate_member(self) -> str:
        account_id = self._get_param("accountId")
        self.inspector2_backend.disassociate_member(account_id)
        return json.dumps({"accountId": account_id})

    def get_member(self) -> str:
        account_id = self._get_param("accountId")
        member = self.inspector2_backend.get_member(account_id)
        return json.dumps({"member": member.to_json()})

    def list_tags_for_resource(self) -> str:
        arn = unquote(self.path.split("/tags/")[-1])
        tags = self.inspector2_backend.list_tags_for_resource(arn)
        return json.dumps({"tags": tags})

    def tag_resource(self) -> str:
        resource_arn = unquote(self.path.split("/tags/")[-1])
        tags = self._get_param("tags")
        self.inspector2_backend.tag_resource(resource_arn=resource_arn, tags=tags)
        return "{}"

    def untag_resource(self) -> str:
        resource_arn = unquote(self.path.split("/tags/")[-1])
        tag_keys = self.querystring.get("tagKeys")
        self.inspector2_backend.untag_resource(resource_arn, tag_keys)  # type: ignore
        return "{}"
