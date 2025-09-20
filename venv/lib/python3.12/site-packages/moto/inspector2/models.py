import json
from typing import Any, Dict, Iterable, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition


class FilterResource(BaseModel):
    def __init__(
        self,
        region: str,
        account_id: str,
        name: str,
        reason: Optional[str],
        action: str,
        description: Optional[str],
        filter_criteria: Dict[str, Any],
        backend: "Inspector2Backend",
    ):
        filter_id = mock_random.get_random_hex(10)
        self.owner_id = account_id
        self.arn = f"arn:{get_partition(region)}:inspector2:{region}:{account_id}:owner/{self.owner_id}/filter/{filter_id}"
        self.name = name
        self.reason = reason
        self.action = action
        self.description = description
        self.filter_criteria = filter_criteria
        self.created_at = unix_time()
        self.backend = backend

    def to_json(self) -> Dict[str, Any]:
        return {
            "action": self.action,
            "arn": self.arn,
            "createdAt": self.created_at,
            "criteria": self.filter_criteria,
            "description": self.description,
            "name": self.name,
            "ownerId": self.owner_id,
            "reason": self.reason,
            "tags": self.backend.list_tags_for_resource(self.arn),
        }


class AccountStatus(BaseModel):
    def __init__(self, account_id: str):
        self.account_id = account_id
        self.ec2 = "DISABLED"
        self.ecr = "DISABLED"
        self._lambda = "DISABLED"
        self.lambda_code = "DISABLED"

    def toggle(self, resource_types: List[str], enable: bool) -> None:
        if "EC2" in resource_types:
            self.ec2 = "ENABLED" if enable else "DISABLED"
        if "ECR" in resource_types:
            self.ecr = "ENABLED" if enable else "DISABLED"
        if "LAMBDA" in resource_types:
            self._lambda = "ENABLED" if enable else "DISABLED"
        if "LAMBDA_CODE" in resource_types or "LAMBDACODE" in resource_types:
            self.lambda_code = "ENABLED" if enable else "DISABLED"

    def to_json(self) -> Dict[str, Any]:
        return {
            "accountId": self.account_id,
            "resourceStatus": {
                "ec2": self.ec2,
                "ecr": self.ecr,
                "lambda": self._lambda,
                "lambdaCode": self.lambda_code,
            },
            "status": self._status(),
        }

    def _status(self) -> str:
        return (
            "ENABLED"
            if "ENABLED" in [self.ec2, self.ecr, self._lambda, self.lambda_code]
            else "DISABLED"
        )

    def to_batch_json(self) -> Dict[str, Any]:
        return {
            "accountId": self.account_id,
            "resourceState": {
                "ec2": {"status": self.ec2},
                "ecr": {"status": self.ecr},
                "lambda": {"status": self._lambda},
                "lambdaCode": {"status": self.lambda_code},
            },
            "state": {"status": self._status()},
        }


class Member(BaseModel):
    def __init__(self, account_id: str, admin_account_id: str):
        self.account_id = account_id
        self.admin_account_id = admin_account_id
        self.status = "ENABLED"
        self.updated_at = unix_time()

    def to_json(self) -> Dict[str, Any]:
        return {
            "accountId": self.account_id,
            "delegatedAdminAccountId": self.admin_account_id,
            "relationshipStatus": self.status,
            "updatedAt": self.updated_at,
        }


class Inspector2Backend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.filters: Dict[str, FilterResource] = dict()
        self.admin_accounts: Dict[str, str] = dict()
        self.account_status: Dict[str, AccountStatus] = dict()
        self.members: Dict[str, Member] = dict()
        self.org_config = {
            "ec2": False,
            "ecr": False,
            "lambda": False,
            "lambdaCode": False,
        }
        self.tagger = TaggingService()
        self.findings_queue: List[Any] = []
        self.findings: Dict[str, Any] = {}

    def create_filter(
        self,
        action: str,
        description: str,
        filter_criteria: Dict[str, Any],
        name: str,
        reason: str,
        tags: Dict[str, str],
    ) -> str:
        _filter = FilterResource(
            region=self.region_name,
            account_id=self.account_id,
            action=action,
            description=description,
            filter_criteria=filter_criteria,
            name=name,
            reason=reason,
            backend=self,
        )
        self.filters[_filter.arn] = _filter
        self.tag_resource(_filter.arn, tags)
        return _filter.arn

    def delete_filter(self, arn: str) -> None:
        self.filters.pop(arn, None)

    def list_filters(self, action: str, arns: List[str]) -> Iterable[FilterResource]:
        """
        Pagination is not yet implemented
        """
        return [
            f
            for f in self.filters.values()
            if (arns and f.arn in arns)
            or (action and f.action == action)
            or (not arns and not action)
        ]

    def list_findings(
        self,
        filter_criteria: List[Dict[str, Any]],
        max_results: str,
        next_token: str,
        sort_criteria: str,
    ) -> List[Dict[str, Any]]:
        """
        This call will always return 0 findings by default.

        You can use a dedicated API to override this, by configuring a queue of expected results.

        A request to `list_findings` will take the first result from that queue, and assign it to the provided arguments. Subsequent calls using the same arguments will return the same result. Other requests using a different SQL-query will take the next result from the queue, or return an empty result if the queue is empty.

        Configure this queue by making an HTTP request to `/moto-api/static/inspector2/findings-results`. An example invocation looks like this:

        .. sourcecode:: python

            findings = {
                "results": [
                    [{
                        "awsAccountId": "111122223333",
                        "codeVulnerabilityDetails": {"cwes": ["a"], "detectorId": ".."},
                    }],
                    # .. other findings as required
                ],
                "account_id": "123456789012",  # This is the default - can be omitted
                "region": "us-east-1",  # This is the default - can be omitted
            }
            resp = requests.post(
                "http://motoapi.amazonaws.com/moto-api/static/inspector2/findings-results",
                json=findings,
            )

            inspector2 = boto3.client("inspector2", region_name="us-east-1")
            findings = inspector2.list_findings()["findings"]

        """
        key = f"{json.dumps(filter_criteria)}--{max_results}--{next_token}--{sort_criteria}"
        if key not in self.findings and self.findings_queue:
            self.findings[key] = self.findings_queue.pop(0)
        if key in self.findings:
            return self.findings[key]
        else:
            return []

    def list_delegated_admin_accounts(self) -> Dict[str, str]:
        return self.admin_accounts

    def enable_delegated_admin_account(self, account_id: str) -> None:
        self.admin_accounts[account_id] = "ENABLED"

    def disable_delegated_admin_account(self, account_id: str) -> None:
        self.admin_accounts[account_id] = "DISABLED"

    def describe_organization_configuration(self) -> Dict[str, Any]:
        return {"autoEnable": self.org_config, "maxAccountLimitReached": False}

    def update_organization_configuration(
        self, auto_enable: Dict[str, bool]
    ) -> Dict[str, Any]:
        self.org_config.update(auto_enable)
        return {"autoEnable": self.org_config}

    def disable(
        self, account_ids: List[str], resource_types: List[str]
    ) -> List[Dict[str, Any]]:
        for acct in account_ids:
            if acct not in self.account_status:
                self.account_status[acct] = AccountStatus(acct)
            self.account_status[acct].toggle(resource_types, enable=False)

        return [
            status.to_json()
            for a_id, status in self.account_status.items()
            if a_id in account_ids
        ]

    def enable(
        self, account_ids: List[str], resource_types: List[str]
    ) -> List[Dict[str, Any]]:
        for acct in account_ids:
            if acct not in self.account_status:
                self.account_status[acct] = AccountStatus(acct)
            self.account_status[acct].toggle(resource_types, enable=True)

        return [
            status.to_json()
            for a_id, status in self.account_status.items()
            if a_id in account_ids
        ]

    def batch_get_account_status(self, account_ids: List[str]) -> List[Dict[str, Any]]:
        return [
            status.to_batch_json()
            for a_id, status in self.account_status.items()
            if a_id in account_ids
        ]

    def list_members(self) -> Iterable[Member]:
        return self.members.values()

    def associate_member(self, account_id: str) -> None:
        self.members[account_id] = Member(
            account_id=account_id, admin_account_id=self.account_id
        )

    def disassociate_member(self, account_id: str) -> None:
        self.members[account_id].status = "DISABLED"

    def get_member(self, account_id: str) -> Member:
        return self.members[account_id]

    def tag_resource(self, resource_arn: str, tags: Dict[str, str]) -> None:
        self.tagger.tag_resource(
            resource_arn, TaggingService.convert_dict_to_tags_input(tags)
        )

    def list_tags_for_resource(self, resource_arn: str) -> Dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_arn)

    def untag_resource(self, arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(arn, tag_keys)


inspector2_backends = BackendDict(Inspector2Backend, "inspector2")
