"""Handles incoming cloudtrail requests, invokes methods, returns responses."""

import json
from typing import Any, Dict

from moto.core.responses import BaseResponse

from .exceptions import InvalidParameterCombinationException
from .models import CloudTrailBackend, cloudtrail_backends


class CloudTrailResponse(BaseResponse):
    """Handler for CloudTrail requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="cloudtrail")

    @property
    def cloudtrail_backend(self) -> CloudTrailBackend:
        """Return backend instance specific for this region."""
        return cloudtrail_backends[self.current_account][self.region]

    def create_trail(self) -> str:
        name = self._get_param("Name")
        bucket_name = self._get_param("S3BucketName")
        is_global = self._get_bool_param("IncludeGlobalServiceEvents", True)
        is_multi_region = self._get_bool_param("IsMultiRegionTrail", False)
        if not is_global and is_multi_region:
            raise InvalidParameterCombinationException(
                "Multi-Region trail must include global service events."
            )
        s3_key_prefix = self._get_param("S3KeyPrefix")
        sns_topic_name = self._get_param("SnsTopicName")
        log_validation = self._get_bool_param("EnableLogFileValidation", False)
        is_org_trail = self._get_bool_param("IsOrganizationTrail", False)
        cw_log_group_arn = self._get_param("CloudWatchLogsLogGroupArn")
        cw_role_arn = self._get_param("CloudWatchLogsRoleArn")
        kms_key_id = self._get_param("KmsKeyId")
        tags_list = self._get_param("TagsList", [])
        trail = self.cloudtrail_backend.create_trail(
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
            tags_list,
        )
        return json.dumps(trail.description())

    def get_trail(self) -> str:
        name = self._get_param("Name")
        trail = self.cloudtrail_backend.get_trail(name)
        return json.dumps({"Trail": trail.description()})

    def get_trail_status(self) -> str:
        name = self._get_param("Name")
        status = self.cloudtrail_backend.get_trail_status(name)
        return json.dumps(status.description())

    def describe_trails(self) -> str:
        include_shadow_trails = self._get_bool_param("includeShadowTrails", True)
        trails = self.cloudtrail_backend.describe_trails(include_shadow_trails)
        return json.dumps(
            {"trailList": [t.description(include_region=True) for t in trails]}
        )

    def list_trails(self) -> str:
        all_trails = self.cloudtrail_backend.list_trails()
        return json.dumps({"Trails": [t.short() for t in all_trails]})

    def start_logging(self) -> str:
        name = self._get_param("Name")
        self.cloudtrail_backend.start_logging(name)
        return json.dumps({})

    def stop_logging(self) -> str:
        name = self._get_param("Name")
        self.cloudtrail_backend.stop_logging(name)
        return json.dumps({})

    def delete_trail(self) -> str:
        name = self._get_param("Name")
        self.cloudtrail_backend.delete_trail(name)
        return json.dumps({})

    def update_trail(self) -> str:
        name = self._get_param("Name")
        s3_bucket_name = self._get_param("S3BucketName")
        s3_key_prefix = self._get_param("S3KeyPrefix")
        sns_topic_name = self._get_param("SnsTopicName")
        include_global_service_events = self._get_param("IncludeGlobalServiceEvents")
        is_multi_region_trail = self._get_param("IsMultiRegionTrail")
        enable_log_file_validation = self._get_param("EnableLogFileValidation")
        is_organization_trail = self._get_param("IsOrganizationTrail")
        cw_log_group_arn = self._get_param("CloudWatchLogsLogGroupArn")
        cw_role_arn = self._get_param("CloudWatchLogsRoleArn")
        kms_key_id = self._get_param("KmsKeyId")
        trail = self.cloudtrail_backend.update_trail(
            name=name,
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
        return json.dumps(trail.description())

    def put_event_selectors(self) -> str:
        params = json.loads(self.body)
        trail_name = params.get("TrailName")
        event_selectors = params.get("EventSelectors")
        advanced_event_selectors = params.get("AdvancedEventSelectors")
        (
            trail_arn,
            event_selectors,
            advanced_event_selectors,
        ) = self.cloudtrail_backend.put_event_selectors(
            trail_name=trail_name,
            event_selectors=event_selectors,
            advanced_event_selectors=advanced_event_selectors,
        )
        return json.dumps(
            dict(
                TrailARN=trail_arn,
                EventSelectors=event_selectors,
                AdvancedEventSelectors=advanced_event_selectors,
            )
        )

    def get_event_selectors(self) -> str:
        params = json.loads(self.body)
        trail_name = params.get("TrailName")
        (
            trail_arn,
            event_selectors,
            advanced_event_selectors,
        ) = self.cloudtrail_backend.get_event_selectors(trail_name=trail_name)
        return json.dumps(
            dict(
                TrailARN=trail_arn,
                EventSelectors=event_selectors,
                AdvancedEventSelectors=advanced_event_selectors,
            )
        )

    def add_tags(self) -> str:
        params = json.loads(self.body)
        resource_id = params.get("ResourceId")
        tags_list = params.get("TagsList")
        self.cloudtrail_backend.add_tags(resource_id=resource_id, tags_list=tags_list)
        return json.dumps(dict())

    def remove_tags(self) -> str:
        resource_id = self._get_param("ResourceId")
        tags_list = self._get_param("TagsList")
        self.cloudtrail_backend.remove_tags(
            resource_id=resource_id, tags_list=tags_list
        )
        return json.dumps(dict())

    def list_tags(self) -> str:
        params = json.loads(self.body)
        resource_id_list = params.get("ResourceIdList")
        resource_tag_list = self.cloudtrail_backend.list_tags(
            resource_id_list=resource_id_list
        )
        return json.dumps(dict(ResourceTagList=resource_tag_list))

    def put_insight_selectors(self) -> str:
        trail_name = self._get_param("TrailName")
        insight_selectors = self._get_param("InsightSelectors")
        trail_arn, insight_selectors = self.cloudtrail_backend.put_insight_selectors(
            trail_name=trail_name, insight_selectors=insight_selectors
        )
        return json.dumps(dict(TrailARN=trail_arn, InsightSelectors=insight_selectors))

    def get_insight_selectors(self) -> str:
        trail_name = self._get_param("TrailName")
        trail_arn, insight_selectors = self.cloudtrail_backend.get_insight_selectors(
            trail_name=trail_name
        )
        resp: Dict[str, Any] = {"TrailARN": trail_arn}
        if insight_selectors:
            resp["InsightSelectors"] = insight_selectors
        return json.dumps(resp)
