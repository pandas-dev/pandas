"""EventBridgeSchedulerBackend class with methods for supported APIs."""

import datetime
from typing import Any, Optional, cast

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time, utcfromtimestamp, utcnow
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import (
    ScheduleExists,
    ScheduleGroupNotFound,
    ScheduleNotFound,
    ValidationException,
)


class Schedule(BaseModel):
    def __init__(
        self,
        region: str,
        account_id: str,
        group_name: str,
        name: str,
        description: Optional[str],
        schedule_expression: str,
        schedule_expression_timezone: Optional[str],
        flexible_time_window: dict[str, Any],
        target: dict[str, Any],
        state: Optional[str],
        kms_key_arn: Optional[str],
        start_date: Optional[str],
        end_date: Optional[str],
        action_after_completion: Optional[str],
    ):
        self.name = name
        self.group_name = group_name
        self.description = description
        self.arn = f"arn:{get_partition(region)}:scheduler:{region}:{account_id}:schedule/{group_name}/{name}"
        self.schedule_expression = schedule_expression
        self.schedule_expression_timezone = schedule_expression_timezone or "UTC"
        self.flexible_time_window = flexible_time_window
        self.target = Schedule.validate_target(target)
        self.state = state or "ENABLED"
        self.kms_key_arn = kms_key_arn
        self.start_date = self._validate_start_date(start_date)
        self.end_date = end_date
        self.action_after_completion = action_after_completion
        self.creation_date = self.last_modified_date = unix_time()

    @staticmethod
    def validate_target(target: dict[str, Any]) -> dict[str, Any]:  # type: ignore[misc]
        if "RetryPolicy" not in target:
            target["RetryPolicy"] = {
                "MaximumEventAgeInSeconds": 86400,
                "MaximumRetryAttempts": 185,
            }
        return target

    def _validate_start_date(self, start_date: Optional[str]) -> Optional[str]:
        # `.count("*")` means a recurrence expression
        # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/scheduler/client/create_schedule.html (StartDate parameter)
        if self.schedule_expression.count("*") and start_date is not None:
            start_date_as_dt = utcfromtimestamp(cast(int, start_date))
            now = utcnow()
            if start_date_as_dt < now:
                diff = now - start_date_as_dt
                rule = datetime.timedelta(minutes=5)
                within_rule = diff <= rule
                if not within_rule:
                    raise ValidationException(
                        message="The StartDate you specify cannot be earlier than 5 minutes ago."
                    )
        return start_date

    def to_dict(self, short: bool = False) -> dict[str, Any]:
        dct: dict[str, Any] = {
            "Arn": self.arn,
            "Name": self.name,
            "GroupName": self.group_name,
            "Description": self.description,
            "ScheduleExpression": self.schedule_expression,
            "ScheduleExpressionTimezone": self.schedule_expression_timezone,
            "FlexibleTimeWindow": self.flexible_time_window,
            "Target": self.target,
            "ActionAfterCompletion": self.action_after_completion,
            "State": self.state,
            "KmsKeyArn": self.kms_key_arn,
            "StartDate": self.start_date,
            "EndDate": self.end_date,
            "CreationDate": self.creation_date,
            "LastModificationDate": self.last_modified_date,
        }
        if short:
            dct["Target"] = {"Arn": dct["Target"]["Arn"]}
        return dct

    def update(
        self,
        description: str,
        end_date: str,
        flexible_time_window: dict[str, Any],
        kms_key_arn: str,
        schedule_expression: str,
        schedule_expression_timezone: str,
        start_date: str,
        state: str,
        target: dict[str, Any],
    ) -> None:
        self.schedule_expression = schedule_expression
        self.schedule_expression_timezone = schedule_expression_timezone or "UTC"
        self.flexible_time_window = flexible_time_window
        self.target = Schedule.validate_target(target)
        self.description = description
        self.state = state
        self.kms_key_arn = kms_key_arn
        self.start_date = start_date
        self.end_date = end_date
        self.last_modified_date = unix_time()


class ScheduleGroup(BaseModel):
    def __init__(self, region: str, account_id: str, name: str):
        self.name = name
        self.arn = f"arn:{get_partition(region)}:scheduler:{region}:{account_id}:schedule-group/{name}"
        self.schedules: dict[str, Schedule] = {}
        self.created_on = None if self.name == "default" else unix_time()
        self.last_modified = None if self.name == "default" else unix_time()

    def add_schedule(self, schedule: Schedule) -> None:
        self.schedules[schedule.name] = schedule

    def get_schedule(self, name: str) -> Schedule:
        if name not in self.schedules:
            raise ScheduleNotFound(name)
        return self.schedules[name]

    def delete_schedule(self, name: str) -> None:
        if name not in self.schedules:
            raise ScheduleNotFound(name)
        self.schedules.pop(name)

    def to_dict(self) -> dict[str, Any]:
        return {
            "Arn": self.arn,
            "CreationDate": self.created_on,
            "LastModificationDate": self.last_modified,
            "Name": self.name,
            "State": "ACTIVE",
        }


class EventBridgeSchedulerBackend(BaseBackend):
    """Implementation of EventBridgeScheduler APIs."""

    PAGINATION_MODEL = {
        "list_schedule_groups": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 50,
            "unique_attribute": "arn",
        },
        "list_schedules": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 50,
            "unique_attribute": "arn",
        },
    }

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.schedules: list[Schedule] = []
        self.schedule_groups = {
            "default": ScheduleGroup(
                region=region_name, account_id=account_id, name="default"
            )
        }
        self.tagger = TaggingService()

    def create_schedule(
        self,
        description: str,
        end_date: str,
        flexible_time_window: dict[str, Any],
        group_name: str,
        kms_key_arn: str,
        name: str,
        schedule_expression: str,
        schedule_expression_timezone: str,
        start_date: str,
        state: str,
        target: dict[str, Any],
        action_after_completion: Optional[str],
    ) -> Schedule:
        """
        The ClientToken parameter is not yet implemented
        """
        group = self.schedule_groups[group_name or "default"]
        if name in group.schedules:
            raise ScheduleExists(name)
        schedule = Schedule(
            region=self.region_name,
            account_id=self.account_id,
            group_name=group.name,
            name=name,
            description=description,
            schedule_expression=schedule_expression,
            schedule_expression_timezone=schedule_expression_timezone,
            flexible_time_window=flexible_time_window,
            target=target,
            state=state,
            kms_key_arn=kms_key_arn,
            start_date=start_date,
            end_date=end_date,
            action_after_completion=action_after_completion,
        )
        group.add_schedule(schedule)
        return schedule

    def get_schedule(self, group_name: Optional[str], name: str) -> Schedule:
        group = self.get_schedule_group(group_name)
        return group.get_schedule(name)

    def delete_schedule(self, group_name: Optional[str], name: str) -> None:
        group = self.get_schedule_group(group_name)
        group.delete_schedule(name)

    def update_schedule(
        self,
        description: str,
        end_date: str,
        flexible_time_window: dict[str, Any],
        group_name: str,
        kms_key_arn: str,
        name: str,
        schedule_expression: str,
        schedule_expression_timezone: str,
        start_date: str,
        state: str,
        target: dict[str, Any],
    ) -> Schedule:
        """
        The ClientToken is not yet implemented
        """
        schedule = self.get_schedule(group_name=group_name, name=name)
        schedule.update(
            description=description,
            end_date=end_date,
            flexible_time_window=flexible_time_window,
            kms_key_arn=kms_key_arn,
            schedule_expression=schedule_expression,
            schedule_expression_timezone=schedule_expression_timezone,
            start_date=start_date,
            state=state,
            target=target,
        )
        return schedule

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_schedules(
        self,
        group_names: Optional[str],
        state: Optional[str],
        name_prefix: Optional[str] = None,
    ) -> list[Schedule]:
        """
        The following parameters are not yet implemented: MaxResults, NamePrefix, NextToken
        """
        results = []
        for group in self.schedule_groups.values():
            if not group_names or group.name in group_names:
                for schedule in group.schedules.values():
                    if not state or schedule.state == state:
                        if not name_prefix or schedule.name.startswith(name_prefix):
                            results.append(schedule)
        # Sort by creation date, newest first
        results.sort(key=lambda x: x.creation_date, reverse=True)
        return results

    def create_schedule_group(
        self, name: str, tags: list[dict[str, str]]
    ) -> ScheduleGroup:
        """
        The ClientToken parameter is not yet implemented
        """
        group = ScheduleGroup(
            region=self.region_name, account_id=self.account_id, name=name
        )
        self.schedule_groups[name] = group
        self.tagger.tag_resource(group.arn, tags)
        return group

    def get_schedule_group(self, group_name: Optional[str]) -> ScheduleGroup:
        if (group_name or "default") not in self.schedule_groups:
            raise ScheduleGroupNotFound(group_name or "default")
        return self.schedule_groups[group_name or "default"]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_schedule_groups(
        self, name_prefix: Optional[str] = None
    ) -> list[ScheduleGroup]:
        results = []
        for group in self.schedule_groups.values():
            if not name_prefix or group.name.startswith(name_prefix):
                results.append(group)
        # Sort by:
        # 1. Default group first
        # 2. Then by creation date, newest first
        results.sort(
            key=lambda x: (
                x.name != "default",  # False (default) comes before True (non-default)
                -1 * (x.created_on or 0),  # Sort by creation date desc, handling None
            )
        )
        return results

    def delete_schedule_group(self, name: Optional[str]) -> None:
        self.schedule_groups.pop(name or "default")

    def list_tags_for_resource(
        self, resource_arn: str
    ) -> dict[str, list[dict[str, str]]]:
        return self.tagger.list_tags_for_resource(resource_arn)

    def tag_resource(self, resource_arn: str, tags: list[dict[str, str]]) -> None:
        self.tagger.tag_resource(resource_arn, tags)

    def untag_resource(self, resource_arn: str, tag_keys: list[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)


scheduler_backends = BackendDict(EventBridgeSchedulerBackend, "scheduler")
