from datetime import datetime
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.utilities.utils import get_partition

from .exceptions import DetectorNotFoundException, FilterNotFoundException


class GuardDutyBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.admin_account_ids: List[str] = []
        self.detectors: Dict[str, Detector] = {}

    def create_detector(
        self,
        enable: bool,
        finding_publishing_frequency: str,
        data_sources: Dict[str, Any],
        tags: Dict[str, str],
    ) -> str:
        if finding_publishing_frequency not in [
            "FIFTEEN_MINUTES",
            "ONE_HOUR",
            "SIX_HOURS",
        ]:
            finding_publishing_frequency = "SIX_HOURS"

        detector = Detector(
            account_id=self.account_id,
            region_name=self.region_name,
            created_at=datetime.now(),
            finding_publish_freq=finding_publishing_frequency,
            enabled=enable,
            datasources=data_sources,
            tags=tags,
        )
        self.detectors[detector.id] = detector
        return detector.id

    def create_filter(
        self,
        detector_id: str,
        name: str,
        action: str,
        description: str,
        finding_criteria: Dict[str, Any],
        rank: int,
    ) -> None:
        detector = self.get_detector(detector_id)
        _filter = Filter(name, action, description, finding_criteria, rank)
        detector.add_filter(_filter)

    def delete_detector(self, detector_id: str) -> None:
        self.detectors.pop(detector_id, None)

    def delete_filter(self, detector_id: str, filter_name: str) -> None:
        detector = self.get_detector(detector_id)
        detector.delete_filter(filter_name)

    def enable_organization_admin_account(self, admin_account_id: str) -> None:
        self.admin_account_ids.append(admin_account_id)

    def list_organization_admin_accounts(self) -> List[str]:
        """
        Pagination is not yet implemented
        """
        return self.admin_account_ids

    def list_detectors(self) -> List[str]:
        """
        The MaxResults and NextToken-parameter have not yet been implemented.
        """
        detectorids = []
        for detector in self.detectors:
            detectorids.append(self.detectors[detector].id)
        return detectorids

    def get_detector(self, detector_id: str) -> "Detector":
        if detector_id not in self.detectors:
            raise DetectorNotFoundException
        return self.detectors[detector_id]

    def get_filter(self, detector_id: str, filter_name: str) -> "Filter":
        detector = self.get_detector(detector_id)
        return detector.get_filter(filter_name)

    def update_detector(
        self,
        detector_id: str,
        enable: bool,
        finding_publishing_frequency: str,
        data_sources: Dict[str, Any],
    ) -> None:
        detector = self.get_detector(detector_id)
        detector.update(enable, finding_publishing_frequency, data_sources)

    def update_filter(
        self,
        detector_id: str,
        filter_name: str,
        action: str,
        description: str,
        finding_criteria: Dict[str, Any],
        rank: int,
    ) -> None:
        detector = self.get_detector(detector_id)
        detector.update_filter(
            filter_name,
            action=action,
            description=description,
            finding_criteria=finding_criteria,
            rank=rank,
        )


class Filter(BaseModel):
    def __init__(
        self,
        name: str,
        action: str,
        description: str,
        finding_criteria: Dict[str, Any],
        rank: int,
    ):
        self.name = name
        self.action = action
        self.description = description
        self.finding_criteria = finding_criteria
        self.rank = rank or 1

    def update(
        self,
        action: Optional[str],
        description: Optional[str],
        finding_criteria: Optional[Dict[str, Any]],
        rank: Optional[int],
    ) -> None:
        if action is not None:
            self.action = action
        if description is not None:
            self.description = description
        if finding_criteria is not None:
            self.finding_criteria = finding_criteria
        if rank is not None:
            self.rank = rank

    def to_json(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "action": self.action,
            "description": self.description,
            "findingCriteria": self.finding_criteria,
            "rank": self.rank,
        }


class Detector(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        created_at: datetime,
        finding_publish_freq: str,
        enabled: bool,
        datasources: Dict[str, Any],
        tags: Dict[str, str],
    ):
        self.id = mock_random.get_random_hex(length=32)
        self.created_at = created_at
        self.finding_publish_freq = finding_publish_freq
        self.service_role = f"arn:{get_partition(region_name)}:iam::{account_id}:role/aws-service-role/guardduty.amazonaws.com/AWSServiceRoleForAmazonGuardDuty"
        self.enabled = enabled
        self.updated_at = created_at
        self.datasources = datasources or {}
        self.tags = tags or {}

        self.filters: Dict[str, Filter] = dict()

    def add_filter(self, _filter: Filter) -> None:
        self.filters[_filter.name] = _filter

    def delete_filter(self, filter_name: str) -> None:
        self.filters.pop(filter_name, None)

    def get_filter(self, filter_name: str) -> Filter:
        if filter_name not in self.filters:
            raise FilterNotFoundException
        return self.filters[filter_name]

    def update_filter(
        self,
        filter_name: str,
        action: str,
        description: str,
        finding_criteria: Dict[str, Any],
        rank: int,
    ) -> None:
        _filter = self.get_filter(filter_name)
        _filter.update(
            action=action,
            description=description,
            finding_criteria=finding_criteria,
            rank=rank,
        )

    def update(
        self,
        enable: bool,
        finding_publishing_frequency: str,
        data_sources: Dict[str, Any],
    ) -> None:
        if enable is not None:
            self.enabled = enable
        if finding_publishing_frequency is not None:
            self.finding_publish_freq = finding_publishing_frequency
        if data_sources is not None:
            self.datasources = data_sources

    def to_json(self) -> Dict[str, Any]:
        data_sources = {
            "cloudTrail": {"status": "DISABLED"},
            "dnsLogs": {"status": "DISABLED"},
            "flowLogs": {"status": "DISABLED"},
            "s3Logs": {
                "status": "ENABLED"
                if (self.datasources.get("s3Logs") or {}).get("enable")
                else "DISABLED"
            },
            "kubernetes": {
                "auditLogs": {
                    "status": "ENABLED"
                    if self.datasources.get("kubernetes", {})
                    .get("auditLogs", {})
                    .get("enable")
                    else "DISABLED"
                }
            },
        }
        return {
            "createdAt": self.created_at.strftime("%Y-%m-%dT%H:%M:%S.%fZ"),
            "findingPublishingFrequency": self.finding_publish_freq,
            "serviceRole": self.service_role,
            "status": "ENABLED" if self.enabled else "DISABLED",
            "updatedAt": self.updated_at.strftime("%Y-%m-%dT%H:%M:%S.%fZ"),
            "dataSources": data_sources,
            "tags": self.tags,
        }


guardduty_backends = BackendDict(GuardDutyBackend, "guardduty")
