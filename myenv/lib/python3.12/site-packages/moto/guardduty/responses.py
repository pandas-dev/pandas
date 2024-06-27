import json
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import GuardDutyBackend, guardduty_backends


class GuardDutyResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="guardduty")

    @property
    def guardduty_backend(self) -> GuardDutyBackend:
        return guardduty_backends[self.current_account][self.region]

    def create_filter(self) -> str:
        detector_id = self.path.split("/")[-2]
        name = self._get_param("name")
        action = self._get_param("action")
        description = self._get_param("description")
        finding_criteria = self._get_param("findingCriteria")
        rank = self._get_param("rank")

        self.guardduty_backend.create_filter(
            detector_id, name, action, description, finding_criteria, rank
        )
        return json.dumps({"name": name})

    def create_detector(self) -> str:
        enable = self._get_param("enable")
        finding_publishing_frequency = self._get_param("findingPublishingFrequency")
        data_sources = self._get_param("dataSources")
        tags = self._get_param("tags")

        detector_id = self.guardduty_backend.create_detector(
            enable, finding_publishing_frequency, data_sources, tags
        )

        return json.dumps(dict(detectorId=detector_id))

    def delete_detector(self) -> str:
        detector_id = self.path.split("/")[-1]

        self.guardduty_backend.delete_detector(detector_id)
        return "{}"

    def delete_filter(self) -> str:
        detector_id = self.path.split("/")[-3]
        filter_name = unquote(self.path.split("/")[-1])

        self.guardduty_backend.delete_filter(detector_id, filter_name)
        return "{}"

    def enable_organization_admin_account(self) -> str:
        admin_account = self._get_param("adminAccountId")
        self.guardduty_backend.enable_organization_admin_account(admin_account)

        return "{}"

    def list_organization_admin_accounts(self) -> str:
        account_ids = self.guardduty_backend.list_organization_admin_accounts()
        accounts = [
            {"adminAccountId": a, "adminStatus": "ENABLED"} for a in account_ids
        ]

        return json.dumps({"adminAccounts": accounts})

    def list_detectors(self) -> str:
        detector_ids = self.guardduty_backend.list_detectors()

        return json.dumps({"detectorIds": detector_ids})

    def get_detector(self) -> str:
        detector_id = self.path.split("/")[-1]

        detector = self.guardduty_backend.get_detector(detector_id)
        return json.dumps(detector.to_json())

    def get_filter(self) -> str:
        detector_id = self.path.split("/")[-3]
        filter_name = unquote(self.path.split("/")[-1])

        _filter = self.guardduty_backend.get_filter(detector_id, filter_name)
        return json.dumps(_filter.to_json())

    def update_detector(self) -> str:
        detector_id = self.path.split("/")[-1]
        enable = self._get_param("enable")
        finding_publishing_frequency = self._get_param("findingPublishingFrequency")
        data_sources = self._get_param("dataSources")

        self.guardduty_backend.update_detector(
            detector_id, enable, finding_publishing_frequency, data_sources
        )
        return "{}"

    def update_filter(self) -> str:
        detector_id = self.path.split("/")[-3]
        filter_name = unquote(self.path.split("/")[-1])
        action = self._get_param("action")
        description = self._get_param("description")
        finding_criteria = self._get_param("findingCriteria")
        rank = self._get_param("rank")

        self.guardduty_backend.update_filter(
            detector_id,
            filter_name,
            action=action,
            description=description,
            finding_criteria=finding_criteria,
            rank=rank,
        )
        return json.dumps({"name": filter_name})
