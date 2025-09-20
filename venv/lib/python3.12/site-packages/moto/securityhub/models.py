"""SecurityHubBackend class with methods for supported APIs."""

import datetime
from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.exceptions import RESTError
from moto.organizations.exceptions import AWSOrganizationsNotInUseException
from moto.organizations.models import organizations_backends
from moto.securityhub.exceptions import InvalidInputException
from moto.utilities.paginator import paginate


class Finding(BaseModel):
    def __init__(self, finding_id: str, finding_data: Dict[str, Any]):
        self.id = finding_id
        self.data = finding_data

    def as_dict(self) -> Dict[str, Any]:
        return self.data


class SecurityHubBackend(BaseBackend):
    """Implementation of SecurityHub APIs."""

    PAGINATION_MODEL = {
        "get_findings": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "Id",
            "fail_on_invalid_token": True,
        }
    }

    _org_configs: Dict[str, Dict[str, Any]] = {}

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.findings: List[Finding] = []
        self.region_name = region_name
        self.org_backend = organizations_backends[self.account_id]["aws"]

    def _get_org_config(self) -> Dict[str, Any]:
        """Get organization config for the current account."""
        try:
            org = self.org_backend.describe_organization()
            org_id = org["Organization"]["Id"]
        except RESTError:
            raise AWSOrganizationsNotInUseException()

        if org_id not in SecurityHubBackend._org_configs:
            SecurityHubBackend._org_configs[org_id] = {
                "admin_account_id": None,
                "auto_enable": False,
                "auto_enable_standards": "DEFAULT",
                "configuration": {
                    "ConfigurationType": "LOCAL",
                    "Status": "ENABLED",
                    "StatusMessage": "",
                },
            }
        return SecurityHubBackend._org_configs[org_id]

    @paginate(pagination_model=PAGINATION_MODEL)
    def get_findings(
        self,
        filters: Optional[Dict[str, Any]] = None,
        sort_criteria: Optional[List[Dict[str, str]]] = None,
        max_results: Optional[int] = None,
        next_token: Optional[str] = None,
    ) -> List[Dict[str, str]]:
        """
        Returns findings based on optional filters and sort criteria.
        """
        if max_results is not None:
            try:
                max_results = int(max_results)
                if max_results < 1 or max_results > 100:
                    raise InvalidInputException(
                        op="GetFindings",
                        msg="MaxResults must be a number between 1 and 100",
                    )
            except ValueError:
                raise InvalidInputException(
                    op="GetFindings", msg="MaxResults must be a number greater than 0"
                )

        findings = self.findings

        # TODO: Apply filters if provided
        # TODO: Apply sort criteria if provided

        return [f.as_dict() for f in findings]

    def batch_import_findings(
        self, findings: List[Dict[str, Any]]
    ) -> Tuple[int, int, List[Dict[str, Any]]]:
        """
        Import findings in batch to SecurityHub.

        Args:
            findings: List of finding dictionaries to import

        Returns:
            Tuple of (failed_count, success_count, failed_findings)
        """
        failed_count = 0
        success_count = 0
        failed_findings = []

        for finding_data in findings:
            try:
                if (
                    not isinstance(finding_data["Resources"], list)
                    or len(finding_data["Resources"]) == 0
                ):
                    raise InvalidInputException(
                        op="BatchImportFindings",
                        msg="Finding must contain at least one resource in the Resources array",
                    )

                finding_id = finding_data["Id"]

                existing_finding = next(
                    (f for f in self.findings if f.id == finding_id), None
                )

                if existing_finding:
                    existing_finding.data.update(finding_data)
                else:
                    new_finding = Finding(finding_id, finding_data)
                    self.findings.append(new_finding)

                success_count += 1

            except Exception as e:
                failed_count += 1
                failed_findings.append(
                    {
                        "Id": finding_data.get("Id", ""),
                        "ErrorCode": "InvalidInput",
                        "ErrorMessage": str(e),
                    }
                )

        return failed_count, success_count, failed_findings

    def enable_organization_admin_account(self, admin_account_id: str) -> None:
        try:
            org = self.org_backend.describe_organization()
            org_id = org["Organization"]["Id"]
        except RESTError:
            raise AWSOrganizationsNotInUseException()

        if self.account_id != org["Organization"]["MasterAccountId"]:
            raise RESTError(
                "AccessDeniedException",
                "The request was rejected because you don't have sufficient permissions "
                "to perform this operation. The security token included in the request "
                "is for an account that isn't authorized to perform this operation.",
            )

        try:
            self.org_backend.get_account_by_id(admin_account_id)
        except RESTError:
            raise RESTError(
                "InvalidInputException",
                f"The request was rejected because the account {admin_account_id} is not "
                f"a member of organization {org_id}.",
            )

        org_config = self._get_org_config()
        org_config["admin_account_id"] = admin_account_id

    def update_organization_configuration(
        self,
        auto_enable: bool,
        auto_enable_standards: Optional[str] = None,
        organization_configuration: Optional[Dict[str, Any]] = None,
    ) -> None:
        try:
            self.org_backend.describe_organization()
        except RESTError:
            raise RESTError(
                "ResourceNotFoundException",
                "The request was rejected because AWS Organizations is not in use or not "
                "configured for this account.",
            )

        org_config = self._get_org_config()
        if not org_config["admin_account_id"]:
            raise RESTError(
                "ResourceNotFoundException",
                "The request was rejected because no administrator account has been designated.",
            )

        if self.account_id != org_config["admin_account_id"]:
            raise RESTError(
                "AccessDeniedException",
                "The request was rejected because you don't have permission to perform "
                "this action. Only the designated administrator account can update the "
                "organization configuration.",
            )

        if organization_configuration:
            config_type = organization_configuration.get("ConfigurationType")
            if config_type not in ["CENTRAL", "LOCAL"]:
                raise RESTError(
                    "InvalidInputException",
                    "The request was rejected because the ConfigurationType value must be "
                    "either CENTRAL or LOCAL.",
                )

            status = organization_configuration.get("Status")
            if status not in ["PENDING", "ENABLED", "FAILED"]:
                raise RESTError(
                    "InvalidInputException",
                    "The request was rejected because the Status value must be one of "
                    "PENDING, ENABLED, or FAILED.",
                )

            if config_type == "CENTRAL":
                if auto_enable:
                    raise RESTError(
                        "ValidationException",
                        "The request was rejected because AutoEnable must be false when "
                        "ConfigurationType is CENTRAL.",
                    )
                if auto_enable_standards != "NONE":
                    raise RESTError(
                        "ValidationException",
                        "The request was rejected because AutoEnableStandards must be NONE "
                        "when ConfigurationType is CENTRAL.",
                    )

            org_config["configuration"] = organization_configuration

        org_config["auto_enable"] = auto_enable

        if auto_enable_standards is not None:
            if auto_enable_standards not in ["NONE", "DEFAULT"]:
                raise RESTError(
                    "InvalidInputException",
                    "The request was rejected because AutoEnableStandards must be either "
                    "NONE or DEFAULT.",
                )
            org_config["auto_enable_standards"] = auto_enable_standards

    def get_administrator_account(self) -> Dict[str, Any]:
        try:
            org = self.org_backend.describe_organization()
            management_account_id = org["Organization"]["MasterAccountId"]
        except RESTError:
            return {}

        org_config = self._get_org_config()
        admin_account_id = org_config["admin_account_id"]

        if not admin_account_id:
            return {}

        if (
            self.account_id == management_account_id
            or self.account_id == admin_account_id
        ):
            return {}

        return {
            "Administrator": {
                "AccountId": admin_account_id,
                "MemberStatus": "ENABLED",
                "InvitationId": f"invitation-{admin_account_id}",
                "InvitedAt": datetime.datetime.now().isoformat(),
            }
        }

    def describe_organization_configuration(self) -> Dict[str, Any]:
        try:
            self.org_backend.describe_organization()
        except RESTError:
            raise RESTError(
                "AccessDeniedException",
                "You do not have sufficient access to perform this action.",
            )

        org_config = self._get_org_config()
        if not org_config["admin_account_id"]:
            raise RESTError(
                "AccessDeniedException",
                "You do not have sufficient access to perform this action.",
            )

        if self.account_id != org_config["admin_account_id"]:
            raise RESTError(
                "AccessDeniedException",
                "You do not have sufficient access to perform this action.",
            )

        return {
            "AutoEnable": org_config["auto_enable"],
            "MemberAccountLimitReached": False,
            "AutoEnableStandards": org_config["auto_enable_standards"],
            "OrganizationConfiguration": org_config["configuration"],
        }


securityhub_backends = BackendDict(SecurityHubBackend, "securityhub")
