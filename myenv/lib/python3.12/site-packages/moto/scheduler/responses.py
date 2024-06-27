"""Handles incoming scheduler requests, invokes methods, returns responses."""

import json
from typing import Any
from urllib.parse import unquote

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import EventBridgeSchedulerBackend, scheduler_backends


class EventBridgeSchedulerResponse(BaseResponse):
    """Handler for EventBridgeScheduler requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="scheduler")

    @property
    def scheduler_backend(self) -> EventBridgeSchedulerBackend:
        """Return backend instance specific for this region."""
        return scheduler_backends[self.current_account][self.region]

    def create_schedule(self) -> str:
        description = self._get_param("Description")
        end_date = self._get_param("EndDate")
        flexible_time_window = self._get_param("FlexibleTimeWindow")
        group_name = self._get_param("GroupName")
        kms_key_arn = self._get_param("KmsKeyArn")
        name = self.uri.split("/")[-1]
        schedule_expression = self._get_param("ScheduleExpression")
        schedule_expression_timezone = self._get_param("ScheduleExpressionTimezone")
        start_date = self._get_param("StartDate")
        state = self._get_param("State")
        target = self._get_param("Target")
        schedule = self.scheduler_backend.create_schedule(
            description=description,
            end_date=end_date,
            flexible_time_window=flexible_time_window,
            group_name=group_name,
            kms_key_arn=kms_key_arn,
            name=name,
            schedule_expression=schedule_expression,
            schedule_expression_timezone=schedule_expression_timezone,
            start_date=start_date,
            state=state,
            target=target,
        )
        return json.dumps(dict(ScheduleArn=schedule.arn))

    def get_schedule(self) -> str:
        group_name = self._get_param("groupName")
        full_url = self.uri.split("?")[0]
        name = full_url.split("/")[-1]
        schedule = self.scheduler_backend.get_schedule(group_name, name)
        return json.dumps(schedule.to_dict())

    def delete_schedule(self) -> str:
        group_name = self._get_param("groupName")
        name = self.uri.split("?")[0].split("/")[-1]
        self.scheduler_backend.delete_schedule(group_name, name)
        return "{}"

    def update_schedule(self) -> str:
        group_name = self._get_param("GroupName")
        name = self.uri.split("?")[0].split("/")[-1]
        description = self._get_param("Description")
        end_date = self._get_param("EndDate")
        flexible_time_window = self._get_param("FlexibleTimeWindow")
        kms_key_arn = self._get_param("KmsKeyArn")
        schedule_expression = self._get_param("ScheduleExpression")
        schedule_expression_timezone = self._get_param("ScheduleExpressionTimezone")
        start_date = self._get_param("StartDate")
        state = self._get_param("State")
        target = self._get_param("Target")
        schedule = self.scheduler_backend.update_schedule(
            description=description,
            end_date=end_date,
            flexible_time_window=flexible_time_window,
            group_name=group_name,
            kms_key_arn=kms_key_arn,
            name=name,
            schedule_expression=schedule_expression,
            schedule_expression_timezone=schedule_expression_timezone,
            start_date=start_date,
            state=state,
            target=target,
        )
        return json.dumps(dict(ScheduleArn=schedule.arn))

    def list_schedules(self) -> str:
        group_names = self.querystring.get("ScheduleGroup")
        state = self._get_param("State")
        schedules = self.scheduler_backend.list_schedules(group_names, state)
        return json.dumps({"Schedules": [sch.to_dict(short=True) for sch in schedules]})

    def create_schedule_group(self) -> str:
        name = self._get_param("Name")
        tags = self._get_param("Tags")
        schedule_group = self.scheduler_backend.create_schedule_group(
            name=name,
            tags=tags,
        )
        return json.dumps(dict(ScheduleGroupArn=schedule_group.arn))

    def get_schedule_group(self) -> str:
        group_name = self.uri.split("?")[0].split("/")[-1]
        group = self.scheduler_backend.get_schedule_group(group_name)
        return json.dumps(group.to_dict())

    def delete_schedule_group(self) -> str:
        group_name = self.uri.split("?")[0].split("/")[-1]
        self.scheduler_backend.delete_schedule_group(group_name)
        return "{}"

    def list_schedule_groups(self) -> str:
        schedule_groups = self.scheduler_backend.list_schedule_groups()
        return json.dumps(dict(ScheduleGroups=[sg.to_dict() for sg in schedule_groups]))

    def list_tags_for_resource(self) -> TYPE_RESPONSE:
        resource_arn = unquote(self.uri.split("/tags/")[-1])
        tags = self.scheduler_backend.list_tags_for_resource(resource_arn)
        return 200, {}, json.dumps(tags)

    def tag_resource(self) -> TYPE_RESPONSE:
        resource_arn = unquote(self.uri.split("/tags/")[-1])
        tags = json.loads(self.body)["Tags"]
        self.scheduler_backend.tag_resource(resource_arn, tags)
        return 200, {}, "{}"

    def untag_resource(self) -> TYPE_RESPONSE:
        resource_arn = unquote(self.uri.split("?")[0].split("/tags/")[-1])
        tag_keys = self.querystring.get("TagKeys")
        self.scheduler_backend.untag_resource(resource_arn, tag_keys)  # type: ignore
        return 200, {}, "{}"

    def tags(self, request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:
        super().setup_class(request, full_url, headers)
        if request.method == "POST":
            return self.tag_resource()
        elif request.method == "DELETE":
            return self.untag_resource()
        else:
            return self.list_tags_for_resource()
