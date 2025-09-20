import re
import time
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_without_milliseconds, utcnow
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import (
    InsufficientSnsTopicPolicyException,
    S3BucketDoesNotExistException,
    TrailNameInvalidChars,
    TrailNameNotEndingCorrectly,
    TrailNameNotStartingCorrectly,
    TrailNameTooLong,
    TrailNameTooShort,
    TrailNotFoundException,
)


def datetime2int(date: datetime) -> int:
    return int(time.mktime(date.timetuple()))


class TrailStatus:
    def __init__(self) -> None:
        self.is_logging = False
        self.latest_delivery_time: Optional[int] = None
        self.latest_delivery_attempt: Optional[str] = ""
        self.started: Optional[datetime] = None
        self.stopped: Optional[datetime] = None

    def start_logging(self) -> None:
        self.is_logging = True
        self.started = utcnow()
        self.latest_delivery_time = datetime2int(utcnow())
        self.latest_delivery_attempt = iso_8601_datetime_without_milliseconds(utcnow())

    def stop_logging(self) -> None:
        self.is_logging = False
        self.stopped = utcnow()

    def description(self) -> Dict[str, Any]:
        if self.is_logging:
            self.latest_delivery_time = datetime2int(utcnow())
            self.latest_delivery_attempt = iso_8601_datetime_without_milliseconds(
                utcnow()
            )
        desc: Dict[str, Any] = {
            "IsLogging": self.is_logging,
            "LatestDeliveryAttemptTime": self.latest_delivery_attempt,
            "LatestNotificationAttemptTime": "",
            "LatestNotificationAttemptSucceeded": "",
            "LatestDeliveryAttemptSucceeded": "",
            "TimeLoggingStarted": "",
            "TimeLoggingStopped": "",
        }
        if self.started:
            desc["StartLoggingTime"] = datetime2int(self.started)
            desc["TimeLoggingStarted"] = iso_8601_datetime_without_milliseconds(
                self.started
            )
            desc["LatestDeliveryTime"] = self.latest_delivery_time
        if self.stopped:
            desc["StopLoggingTime"] = datetime2int(self.stopped)
            desc["TimeLoggingStopped"] = iso_8601_datetime_without_milliseconds(
                self.stopped
            )
        return desc


class Trail(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        trail_name: str,
        bucket_name: str,
        s3_key_prefix: str,
        sns_topic_name: str,
        is_global: bool,
        is_multi_region: bool,
        log_validation: bool,
        is_org_trail: bool,
        cw_log_group_arn: str,
        cw_role_arn: str,
        kms_key_id: str,
    ):
        self.account_id = account_id
        self.region_name = region_name
        self.partition = get_partition(region_name)
        self.trail_name = trail_name
        self.bucket_name = bucket_name
        self.s3_key_prefix = s3_key_prefix
        self.sns_topic_name = sns_topic_name
        self.is_multi_region = is_multi_region
        self.log_validation = log_validation
        self.is_org_trail = is_org_trail
        self.include_global_service_events = is_global
        self.cw_log_group_arn = cw_log_group_arn
        self.cw_role_arn = cw_role_arn
        self.kms_key_id = kms_key_id
        self.check_name()
        self.check_bucket_exists()
        self.check_topic_exists()
        self.status = TrailStatus()
        self.event_selectors: List[Dict[str, Any]] = list()
        self.advanced_event_selectors: List[Dict[str, Any]] = list()
        self.insight_selectors: List[Dict[str, str]] = list()

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.region_name)}:cloudtrail:{self.region_name}:{self.account_id}:trail/{self.trail_name}"

    @property
    def topic_arn(self) -> Optional[str]:
        if self.sns_topic_name:
            return f"arn:{get_partition(self.region_name)}:sns:{self.region_name}:{self.account_id}:{self.sns_topic_name}"
        return None

    def check_name(self) -> None:
        if len(self.trail_name) < 3:
            raise TrailNameTooShort(actual_length=len(self.trail_name))
        if len(self.trail_name) > 128:
            raise TrailNameTooLong(actual_length=len(self.trail_name))
        if not re.match("^[0-9a-zA-Z]{1}.+$", self.trail_name):
            raise TrailNameNotStartingCorrectly()
        if not re.match(r".+[0-9a-zA-Z]{1}$", self.trail_name):
            raise TrailNameNotEndingCorrectly()
        if not re.match(r"^[.\-_0-9a-zA-Z]+$", self.trail_name):
            raise TrailNameInvalidChars()

    def check_bucket_exists(self) -> None:
        from moto.s3.models import s3_backends

        try:
            s3_backends[self.account_id][self.partition].get_bucket(self.bucket_name)
        except Exception:
            raise S3BucketDoesNotExistException(
                f"S3 bucket {self.bucket_name} does not exist!"
            )

    def check_topic_exists(self) -> None:
        if self.topic_arn:
            from moto.sns import sns_backends

            sns_backend = sns_backends[self.account_id][self.region_name]
            try:
                sns_backend.get_topic(self.topic_arn)
            except Exception:
                raise InsufficientSnsTopicPolicyException(
                    "SNS Topic does not exist or the topic policy is incorrect!"
                )

    def start_logging(self) -> None:
        self.status.start_logging()

    def stop_logging(self) -> None:
        self.status.stop_logging()

    def put_event_selectors(
        self,
        event_selectors: List[Dict[str, Any]],
        advanced_event_selectors: List[Dict[str, Any]],
    ) -> None:
        if event_selectors:
            self.event_selectors = event_selectors
        elif advanced_event_selectors:
            self.event_selectors = []
            self.advanced_event_selectors = advanced_event_selectors

    def get_event_selectors(self) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        return self.event_selectors, self.advanced_event_selectors

    def put_insight_selectors(self, insight_selectors: List[Dict[str, str]]) -> None:
        self.insight_selectors.extend(insight_selectors)

    def get_insight_selectors(self) -> List[Dict[str, str]]:
        return self.insight_selectors

    def update(
        self,
        s3_bucket_name: Optional[str],
        s3_key_prefix: Optional[str],
        sns_topic_name: Optional[str],
        include_global_service_events: Optional[bool],
        is_multi_region_trail: Optional[bool],
        enable_log_file_validation: Optional[bool],
        is_organization_trail: Optional[bool],
        cw_log_group_arn: Optional[str],
        cw_role_arn: Optional[str],
        kms_key_id: Optional[str],
    ) -> None:
        if s3_bucket_name is not None:
            self.bucket_name = s3_bucket_name
        if s3_key_prefix is not None:
            self.s3_key_prefix = s3_key_prefix
        if sns_topic_name is not None:
            self.sns_topic_name = sns_topic_name
        if include_global_service_events is not None:
            self.include_global_service_events = include_global_service_events
        if is_multi_region_trail is not None:
            self.is_multi_region = is_multi_region_trail
        if enable_log_file_validation is not None:
            self.log_validation = enable_log_file_validation
        if is_organization_trail is not None:
            self.is_org_trail = is_organization_trail
        if cw_log_group_arn is not None:
            self.cw_log_group_arn = cw_log_group_arn
        if cw_role_arn is not None:
            self.cw_role_arn = cw_role_arn
        if kms_key_id is not None:
            self.kms_key_id = kms_key_id

    def short(self) -> Dict[str, str]:
        return {
            "Name": self.trail_name,
            "TrailARN": self.arn,
            "HomeRegion": self.region_name,
        }

    def description(self, include_region: bool = False) -> Dict[str, Any]:
        desc = {
            "Name": self.trail_name,
            "S3BucketName": self.bucket_name,
            "IncludeGlobalServiceEvents": self.include_global_service_events,
            "IsMultiRegionTrail": self.is_multi_region,
            "TrailARN": self.arn,
            "LogFileValidationEnabled": self.log_validation,
            "IsOrganizationTrail": self.is_org_trail,
            "HasCustomEventSelectors": False,
            "HasInsightSelectors": False,
            "CloudWatchLogsLogGroupArn": self.cw_log_group_arn,
            "CloudWatchLogsRoleArn": self.cw_role_arn,
            "KmsKeyId": self.kms_key_id,
        }
        if self.s3_key_prefix is not None:
            desc["S3KeyPrefix"] = self.s3_key_prefix
        if self.sns_topic_name is not None:
            desc["SnsTopicName"] = self.sns_topic_name
            desc["SnsTopicARN"] = self.topic_arn
        if include_region:
            desc["HomeRegion"] = self.region_name
        return desc


class CloudTrailBackend(BaseBackend):
    """Implementation of CloudTrail APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.trails: Dict[str, Trail] = dict()
        self.tagging_service = TaggingService(tag_name="TagsList")

    def create_trail(
        self,
        name: str,
        bucket_name: str,
        s3_key_prefix: str,
        sns_topic_name: str,
        is_global: bool,
        is_multi_region: bool,
        log_validation: bool,
        is_org_trail: bool,
        cw_log_group_arn: str,
        cw_role_arn: str,
        kms_key_id: str,
        tags_list: List[Dict[str, str]],
    ) -> Trail:
        trail = Trail(
            self.account_id,
            self.region_name,
            name,
            bucket_name,
            s3_key_prefix,
            sns_topic_name,
            is_global,
            is_multi_region,
            log_validation,
            is_org_trail,
            cw_log_group_arn,
            cw_role_arn,
            kms_key_id,
        )
        self.trails[name] = trail
        self.tagging_service.tag_resource(trail.arn, tags_list)
        return trail

    def get_trail(self, name_or_arn: str) -> Trail:
        if len(name_or_arn) < 3:
            raise TrailNameTooShort(actual_length=len(name_or_arn))
        if name_or_arn in self.trails:
            return self.trails[name_or_arn]
        for trail in self.trails.values():
            if trail.arn == name_or_arn:
                return trail
        raise TrailNotFoundException(account_id=self.account_id, name=name_or_arn)

    def get_trail_status(self, name: str) -> TrailStatus:
        if len(name) < 3:
            raise TrailNameTooShort(actual_length=len(name))

        all_trails = self.describe_trails(include_shadow_trails=True)
        trail = next(
            (
                trail
                for trail in all_trails
                if trail.trail_name == name or trail.arn == name
            ),
            None,
        )
        if not trail:
            # This particular method returns the ARN as part of the error message
            arn = f"arn:{get_partition(self.region_name)}:cloudtrail:{self.region_name}:{self.account_id}:trail/{name}"
            raise TrailNotFoundException(account_id=self.account_id, name=arn)
        return trail.status

    def describe_trails(self, include_shadow_trails: bool) -> Iterable[Trail]:
        all_trails = []
        if include_shadow_trails:
            current_account = cloudtrail_backends[self.account_id]
            for backend in current_account.values():
                for trail in backend.trails.values():
                    if trail.is_multi_region or trail.region_name == self.region_name:
                        all_trails.append(trail)
        else:
            all_trails.extend(self.trails.values())
        return all_trails

    def list_trails(self) -> Iterable[Trail]:
        return self.describe_trails(include_shadow_trails=True)

    def start_logging(self, name: str) -> None:
        trail = self.trails[name]
        trail.start_logging()

    def stop_logging(self, name: str) -> None:
        trail = self.trails[name]
        trail.stop_logging()

    def delete_trail(self, name: str) -> None:
        if name in self.trails:
            del self.trails[name]

    def update_trail(
        self,
        name: str,
        s3_bucket_name: str,
        s3_key_prefix: str,
        sns_topic_name: str,
        include_global_service_events: bool,
        is_multi_region_trail: bool,
        enable_log_file_validation: bool,
        is_organization_trail: bool,
        cw_log_group_arn: str,
        cw_role_arn: str,
        kms_key_id: str,
    ) -> Trail:
        trail = self.get_trail(name_or_arn=name)
        trail.update(
            s3_bucket_name=s3_bucket_name,
            s3_key_prefix=s3_key_prefix,
            sns_topic_name=sns_topic_name,
            include_global_service_events=include_global_service_events,
            is_multi_region_trail=is_multi_region_trail,
            enable_log_file_validation=enable_log_file_validation,
            is_organization_trail=is_organization_trail,
            cw_log_group_arn=cw_log_group_arn,
            cw_role_arn=cw_role_arn,
            kms_key_id=kms_key_id,
        )
        return trail

    def put_event_selectors(
        self,
        trail_name: str,
        event_selectors: List[Dict[str, Any]],
        advanced_event_selectors: List[Dict[str, Any]],
    ) -> Tuple[str, List[Dict[str, Any]], List[Dict[str, Any]]]:
        trail = self.get_trail(trail_name)
        trail.put_event_selectors(event_selectors, advanced_event_selectors)
        trail_arn = trail.arn
        return trail_arn, event_selectors, advanced_event_selectors

    def get_event_selectors(
        self, trail_name: str
    ) -> Tuple[str, List[Dict[str, Any]], List[Dict[str, Any]]]:
        trail = self.get_trail(trail_name)
        event_selectors, advanced_event_selectors = trail.get_event_selectors()
        return trail.arn, event_selectors, advanced_event_selectors

    def add_tags(self, resource_id: str, tags_list: List[Dict[str, str]]) -> None:
        self.tagging_service.tag_resource(resource_id, tags_list)

    def remove_tags(self, resource_id: str, tags_list: List[Dict[str, str]]) -> None:
        self.tagging_service.untag_resource_using_tags(resource_id, tags_list)

    def list_tags(self, resource_id_list: List[str]) -> List[Dict[str, Any]]:
        """
        Pagination is not yet implemented
        """
        resp: List[Dict[str, Any]] = [{"ResourceId": r_id} for r_id in resource_id_list]
        for item in resp:
            item["TagsList"] = self.tagging_service.list_tags_for_resource(
                item["ResourceId"]
            )["TagsList"]
        return resp

    def put_insight_selectors(
        self, trail_name: str, insight_selectors: List[Dict[str, str]]
    ) -> Tuple[str, List[Dict[str, str]]]:
        trail = self.get_trail(trail_name)
        trail.put_insight_selectors(insight_selectors)
        return trail.arn, insight_selectors

    def get_insight_selectors(
        self, trail_name: str
    ) -> Tuple[str, List[Dict[str, str]]]:
        trail = self.get_trail(trail_name)
        return trail.arn, trail.get_insight_selectors()


cloudtrail_backends = BackendDict(CloudTrailBackend, "cloudtrail")
