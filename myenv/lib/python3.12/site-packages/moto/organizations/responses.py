import json
from typing import Any, Dict

from moto.core.responses import BaseResponse

from .models import OrganizationsBackend, organizations_backends


class OrganizationsResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="organizations")

    @property
    def organizations_backend(self) -> OrganizationsBackend:
        return organizations_backends[self.current_account][self.partition]

    @property
    def request_params(self) -> Dict[str, Any]:  # type: ignore[misc]
        try:
            return json.loads(self.body)
        except ValueError:
            return {}

    def _get_param(self, param_name: str, if_none: Any = None) -> Any:
        return self.request_params.get(param_name, if_none)

    def create_organization(self) -> str:
        return json.dumps(
            self.organizations_backend.create_organization(
                region=self.region, **self.request_params
            )
        )

    def describe_organization(self) -> str:
        return json.dumps(self.organizations_backend.describe_organization())

    def delete_organization(self) -> str:
        self.organizations_backend.delete_organization()
        return "{}"

    def list_roots(self) -> str:
        return json.dumps(self.organizations_backend.list_roots())

    def create_organizational_unit(self) -> str:
        return json.dumps(
            self.organizations_backend.create_organizational_unit(**self.request_params)
        )

    def delete_organizational_unit(self) -> str:
        self.organizations_backend.delete_organizational_unit(**self.request_params)
        return "{}"

    def update_organizational_unit(self) -> str:
        return json.dumps(
            self.organizations_backend.update_organizational_unit(**self.request_params)
        )

    def describe_organizational_unit(self) -> str:
        return json.dumps(
            self.organizations_backend.describe_organizational_unit(
                **self.request_params
            )
        )

    def list_organizational_units_for_parent(self) -> str:
        max_results = self._get_int_param("MaxResults")
        next_token = self._get_param("NextToken")
        parent_id = self._get_param("ParentId")
        (
            ous,
            next_token,
        ) = self.organizations_backend.list_organizational_units_for_parent(
            max_results=max_results, next_token=next_token, parent_id=parent_id
        )
        response = {"OrganizationalUnits": ous}
        if next_token:
            response["NextToken"] = next_token
        return json.dumps(response)

    def list_parents(self) -> str:
        return json.dumps(
            self.organizations_backend.list_parents(**self.request_params)
        )

    def create_account(self) -> str:
        return json.dumps(
            self.organizations_backend.create_account(**self.request_params)
        )

    def close_account(self) -> str:
        self.organizations_backend.close_account(**self.request_params)
        return "{}"

    def describe_account(self) -> str:
        return json.dumps(
            self.organizations_backend.describe_account(**self.request_params)
        )

    def describe_create_account_status(self) -> str:
        return json.dumps(
            self.organizations_backend.describe_create_account_status(
                **self.request_params
            )
        )

    def list_create_account_status(self) -> str:
        return json.dumps(
            self.organizations_backend.list_create_account_status(**self.request_params)
        )

    def list_accounts(self) -> str:
        max_results = self._get_int_param("MaxResults")
        next_token = self._get_param("NextToken")
        accounts, next_token = self.organizations_backend.list_accounts(
            max_results=max_results, next_token=next_token
        )
        response = {"Accounts": accounts}
        if next_token:
            response["NextToken"] = next_token
        return json.dumps(response)

    def list_accounts_for_parent(self) -> str:
        max_results = self._get_int_param("MaxResults")
        next_token = self._get_param("NextToken")
        parent_id = self._get_param("ParentId")
        accounts, next_token = self.organizations_backend.list_accounts_for_parent(
            max_results=max_results, next_token=next_token, parent_id=parent_id
        )
        response = {"Accounts": accounts}
        if next_token:
            response["NextToken"] = next_token
        return json.dumps(response)

    def move_account(self) -> str:
        self.organizations_backend.move_account(**self.request_params)
        return "{}"

    def list_children(self) -> str:
        return json.dumps(
            self.organizations_backend.list_children(**self.request_params)
        )

    def create_policy(self) -> str:
        return json.dumps(
            self.organizations_backend.create_policy(**self.request_params)
        )

    def describe_policy(self) -> str:
        return json.dumps(
            self.organizations_backend.describe_policy(**self.request_params)
        )

    def update_policy(self) -> str:
        return json.dumps(
            self.organizations_backend.update_policy(**self.request_params)
        )

    def attach_policy(self) -> str:
        self.organizations_backend.attach_policy(**self.request_params)
        return "{}"

    def list_policies(self) -> str:
        return json.dumps(self.organizations_backend.list_policies())

    def delete_policy(self) -> str:
        self.organizations_backend.delete_policy(**self.request_params)
        return json.dumps({})

    def list_policies_for_target(self) -> str:
        return json.dumps(
            self.organizations_backend.list_policies_for_target(**self.request_params)
        )

    def list_targets_for_policy(self) -> str:
        return json.dumps(
            self.organizations_backend.list_targets_for_policy(**self.request_params)
        )

    def tag_resource(self) -> str:
        self.organizations_backend.tag_resource(**self.request_params)
        return "{}"

    def list_tags_for_resource(self) -> str:
        return json.dumps(
            self.organizations_backend.list_tags_for_resource(**self.request_params)
        )

    def untag_resource(self) -> str:
        self.organizations_backend.untag_resource(**self.request_params)
        return "{}"

    def enable_aws_service_access(self) -> str:
        self.organizations_backend.enable_aws_service_access(**self.request_params)
        return "{}"

    def list_aws_service_access_for_organization(self) -> str:
        return json.dumps(
            self.organizations_backend.list_aws_service_access_for_organization()
        )

    def disable_aws_service_access(self) -> str:
        self.organizations_backend.disable_aws_service_access(**self.request_params)
        return "{}"

    def register_delegated_administrator(self) -> str:
        self.organizations_backend.register_delegated_administrator(
            **self.request_params
        )
        return "{}"

    def list_delegated_administrators(self) -> str:
        return json.dumps(
            self.organizations_backend.list_delegated_administrators(
                **self.request_params
            )
        )

    def list_delegated_services_for_account(self) -> str:
        return json.dumps(
            self.organizations_backend.list_delegated_services_for_account(
                **self.request_params
            )
        )

    def deregister_delegated_administrator(self) -> str:
        self.organizations_backend.deregister_delegated_administrator(
            **self.request_params
        )
        return "{}"

    def enable_policy_type(self) -> str:
        return json.dumps(
            self.organizations_backend.enable_policy_type(**self.request_params)
        )

    def disable_policy_type(self) -> str:
        return json.dumps(
            self.organizations_backend.disable_policy_type(**self.request_params)
        )

    def detach_policy(self) -> str:
        self.organizations_backend.detach_policy(**self.request_params)
        return "{}"

    def remove_account_from_organization(self) -> str:
        self.organizations_backend.remove_account_from_organization(
            **self.request_params
        )
        return "{}"
