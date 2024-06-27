import datetime
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.moto_api._internal import mock_random as random
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.utilities.utils import load_resource

checks_json = "resources/describe_trusted_advisor_checks.json"
ADVISOR_CHECKS = load_resource(__name__, checks_json)


class SupportCase(ManagedState):
    def __init__(self, **kwargs: Any):
        # Configure ManagedState
        super().__init__(
            "support::case",
            transitions=[
                ("opened", "pending-customer-action"),
                ("pending-customer-action", "reopened"),
                ("reopened", "resolved"),
                ("resolved", "unassigned"),
                ("unassigned", "work-in-progress"),
                ("work-in-progress", "opened"),
            ],
        )
        self.case_id = kwargs.get("case_id")
        self.display_id = "foo_display_id"
        self.subject = kwargs.get("subject")
        self.service_code = kwargs.get("service_code")
        self.category_code = kwargs.get("category_code")
        self.severity_code = kwargs.get("severity_code")
        self.submitted_by = "moto@moto.com"
        self.time_created = self.get_datetime()
        self.attachment_set_id = kwargs.get("attachment_set_id")
        self.communication_body = kwargs.get("communication_body")
        self.language = kwargs.get("language")
        self.cc_email_addresses = kwargs.get("cc_email_addresses")
        self.communications = {
            "recentCommunications": {
                "communications": [
                    {
                        "caseId": self.case_id,
                        "body": self.communication_body,
                        "submittedBy": self.submitted_by,
                        "timeCreated": self.get_datetime(),
                        "attachmentSet": [
                            {
                                "attachmentId": self.attachment_set_id,
                                "fileName": "support_file.txt",
                            },
                        ],
                    }
                ],
                "nextToken": "foo_next_token",
            }
        }

    def get_datetime(self) -> str:
        return str(datetime.datetime.now().isoformat())


class SupportBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.check_status: Dict[str, str] = {}
        self.cases: Dict[str, SupportCase] = {}

    def describe_trusted_advisor_checks(self) -> List[Dict[str, Any]]:
        """
        The Language-parameter is not yet implemented
        """
        # The checks are a static response
        checks = ADVISOR_CHECKS["checks"]
        return checks

    def refresh_trusted_advisor_check(self, check_id: str) -> Dict[str, Any]:
        self.advance_check_status(check_id)
        return {
            "status": {
                "checkId": check_id,
                "status": self.check_status[check_id],
                "millisUntilNextRefreshable": 123,
            }
        }

    def advance_check_status(self, check_id: str) -> None:
        """
        Fake an advancement through statuses on refreshing TA checks
        """
        if check_id not in self.check_status:
            self.check_status[check_id] = "none"

        elif self.check_status[check_id] == "none":
            self.check_status[check_id] = "enqueued"

        elif self.check_status[check_id] == "enqueued":
            self.check_status[check_id] = "processing"

        elif self.check_status[check_id] == "processing":
            self.check_status[check_id] = "success"

        elif self.check_status[check_id] == "success":
            self.check_status[check_id] = "abandoned"

        elif self.check_status[check_id] == "abandoned":
            self.check_status[check_id] = "none"

    def advance_case_status(self, case_id: str) -> None:
        """
        Fake an advancement through case statuses
        """
        self.cases[case_id].advance()

    def advance_case_severity_codes(self, case_id: str) -> None:
        """
        Fake an advancement through case status severities
        """
        if self.cases[case_id].severity_code == "low":
            self.cases[case_id].severity_code = "normal"

        elif self.cases[case_id].severity_code == "normal":
            self.cases[case_id].severity_code = "high"

        elif self.cases[case_id].severity_code == "high":
            self.cases[case_id].severity_code = "urgent"

        elif self.cases[case_id].severity_code == "urgent":
            self.cases[case_id].severity_code = "critical"

        elif self.cases[case_id].severity_code == "critical":
            self.cases[case_id].severity_code = "low"

    def resolve_case(self, case_id: str) -> Dict[str, Optional[str]]:
        self.advance_case_status(case_id)

        return {
            "initialCaseStatus": self.cases[case_id].status,
            "finalCaseStatus": "resolved",
        }

    # persist case details to self.cases
    def create_case(
        self,
        subject: str,
        service_code: str,
        severity_code: str,
        category_code: str,
        communication_body: str,
        cc_email_addresses: List[str],
        language: str,
        attachment_set_id: str,
    ) -> Dict[str, str]:
        """
        The IssueType-parameter is not yet implemented
        """
        # Random case ID
        random_case_id = "".join(
            random.choice("0123456789ABCDEFGHIJKLMabcdefghijklm") for i in range(16)
        )
        case_id = f"case-12345678910-2020-{random_case_id}"
        case = SupportCase(
            case_id=case_id,
            subject=subject,
            service_code=service_code,
            severity_code=severity_code,
            category_code=category_code,
            communication_body=communication_body,
            cc_email_addresses=cc_email_addresses,
            language=language,
            attachment_set_id=attachment_set_id,
        )
        self.cases[case_id] = case

        return {"caseId": case_id}

    def describe_cases(
        self,
        case_id_list: List[str],
        include_resolved_cases: bool,
        next_token: Optional[str],
        include_communications: bool,
    ) -> Dict[str, Any]:
        """
        The following parameters have not yet been implemented:
        DisplayID, AfterTime, BeforeTime, MaxResults, Language
        """
        cases = []
        requested_case_ids = case_id_list or self.cases.keys()

        for case in requested_case_ids:
            self.advance_case_status(case)
            self.advance_case_severity_codes(case)
            formatted_case = {
                "caseId": self.cases[case].case_id,
                "displayId": self.cases[case].display_id,
                "subject": self.cases[case].subject,
                "status": self.cases[case].status,
                "serviceCode": self.cases[case].service_code,
                "categoryCode": self.cases[case].category_code,
                "severityCode": self.cases[case].severity_code,
                "submittedBy": self.cases[case].submitted_by,
                "timeCreated": self.cases[case].time_created,
                "ccEmailAddresses": self.cases[case].cc_email_addresses,
                "language": self.cases[case].language,
            }

            if include_communications:
                formatted_case.update(self.cases[case].communications)

            if (
                include_resolved_cases is False
                and formatted_case["status"] == "resolved"
            ):
                continue

            cases.append(formatted_case)
        return {"cases": cases, "nextToken": next_token}


support_backends = BackendDict(
    SupportBackend, "support", use_boto3_regions=False, additional_regions=["us-east-1"]
)
