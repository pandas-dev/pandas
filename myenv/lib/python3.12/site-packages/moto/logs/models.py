from datetime import datetime, timedelta
from gzip import compress as gzip_compress
from typing import Any, Dict, Iterable, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import unix_time_millis, utcnow
from moto.logs.exceptions import (
    InvalidParameterException,
    LimitExceededException,
    ResourceAlreadyExistsException,
    ResourceNotFoundException,
)
from moto.logs.logs_query import execute_query
from moto.logs.metric_filters import MetricFilters
from moto.moto_api._internal import mock_random
from moto.s3.models import MissingBucket, s3_backends
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .utils import PAGINATION_MODEL, EventMessageFilter

MAX_RESOURCE_POLICIES_PER_REGION = 10


class Destination(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        destination_name: str,
        role_arn: str,
        target_arn: str,
        access_policy: Optional[str] = None,
    ):
        self.access_policy = access_policy
        self.arn = f"arn:{get_partition(region)}:logs:{region}:{account_id}:destination:{destination_name}"
        self.creation_time = int(unix_time_millis())
        self.destination_name = destination_name
        self.role_arn = role_arn
        self.target_arn = target_arn

    def to_dict(self) -> Dict[str, Any]:
        return {
            "accessPolicy": self.access_policy,
            "arn": self.arn,
            "creationTime": self.creation_time,
            "destinationName": self.destination_name,
            "roleArn": self.role_arn,
            "targetArn": self.target_arn,
        }


class LogQuery(BaseModel):
    def __init__(
        self,
        query_id: str,
        start_time: int,
        end_time: int,
        query: str,
        log_groups: List["LogGroup"],
    ):
        self.query_id = query_id
        self.start_time = start_time
        self.end_time = end_time
        self.query = query
        self.log_group_names = [lg.name for lg in log_groups]
        self.create_time = unix_time_millis()
        self.status = "Running"
        self.results = execute_query(
            log_groups=log_groups, query=query, start_time=start_time, end_time=end_time
        )
        self.status = "Complete"

    def to_json(self, log_group_name: str) -> Dict[str, Any]:
        return {
            "queryId": self.query_id,
            "queryString": self.query,
            "status": self.status,
            "createTime": self.create_time,
            "logGroupName": log_group_name,
        }

    def to_result_json(self) -> Dict[str, Any]:
        return {
            "results": [
                [{"field": key, "value": val} for key, val in result.items()]
                for result in self.results
            ],
            "status": self.status,
        }


class LogEvent(BaseModel):
    _event_id = 0

    def __init__(self, ingestion_time: int, log_event: Dict[str, Any]):
        self.ingestion_time = ingestion_time
        self.timestamp = log_event["timestamp"]
        self.message = log_event["message"]
        self.event_id = self.__class__._event_id
        self.__class__._event_id += 1
        ""

    def to_filter_dict(self) -> Dict[str, Any]:
        return {
            "eventId": str(self.event_id),
            "ingestionTime": self.ingestion_time,
            # "logStreamName":
            "message": self.message,
            "timestamp": self.timestamp,
        }

    def to_response_dict(self) -> Dict[str, Any]:
        return {
            "ingestionTime": self.ingestion_time,
            "message": self.message,
            "timestamp": self.timestamp,
        }


class LogStream(BaseModel):
    _log_ids = 0

    def __init__(self, log_group: "LogGroup", name: str):
        self.account_id = log_group.account_id
        self.region = log_group.region
        self.log_group = log_group
        self.arn = f"arn:{get_partition(self.region)}:logs:{self.region}:{self.account_id}:log-group:{log_group.name}:log-stream:{name}"
        self.creation_time = int(unix_time_millis())
        self.first_event_timestamp = None
        self.last_event_timestamp = None
        self.last_ingestion_time: Optional[int] = None
        self.log_stream_name = name
        self.stored_bytes = 0
        # I'm  guessing this is token needed for sequenceToken by put_events
        self.upload_sequence_token = 0
        self.events: List[LogEvent] = []

        self.__class__._log_ids += 1

    def _update(self) -> None:
        # events can be empty when stream is described soon after creation
        self.first_event_timestamp = (
            min([x.timestamp for x in self.events]) if self.events else None
        )
        self.last_event_timestamp = (
            max([x.timestamp for x in self.events]) if self.events else None
        )

    def to_describe_dict(self) -> Dict[str, Any]:
        # Compute start and end times
        self._update()

        res = {
            "arn": self.arn,
            "creationTime": self.creation_time,
            "logStreamName": self.log_stream_name,
            "storedBytes": self.stored_bytes,
        }
        if self.events:
            rest = {
                "firstEventTimestamp": self.first_event_timestamp,
                "lastEventTimestamp": self.last_event_timestamp,
                "lastIngestionTime": self.last_ingestion_time,
                "uploadSequenceToken": str(self.upload_sequence_token),
            }
            res.update(rest)
        return res

    def put_log_events(self, log_events: List[Dict[str, Any]]) -> str:
        # TODO: ensure sequence_token
        # TODO: to be thread safe this would need a lock
        self.last_ingestion_time = int(unix_time_millis())
        # TODO: make this match AWS if possible
        self.stored_bytes += sum(
            [len(log_event["message"]) for log_event in log_events]
        )
        events = [
            LogEvent(self.last_ingestion_time, log_event) for log_event in log_events
        ]
        self.events += events
        self.upload_sequence_token += 1

        for subscription_filter in self.log_group.subscription_filters.values():
            service = subscription_filter.destination_arn.split(":")[2]
            formatted_log_events = [
                {
                    "id": event.event_id,
                    "timestamp": event.timestamp,
                    "message": event.message,
                }
                for event in events
            ]
            self._send_log_events(
                service=service,
                destination_arn=subscription_filter.destination_arn,
                filter_name=subscription_filter.name,
                log_events=formatted_log_events,
            )
        return f"{self.upload_sequence_token:056d}"

    def _send_log_events(
        self,
        service: str,
        destination_arn: str,
        filter_name: str,
        log_events: List[Dict[str, Any]],
    ) -> None:
        if service == "lambda":
            from moto.awslambda.utils import get_backend

            get_backend(self.account_id, self.region).send_log_event(
                destination_arn,
                filter_name,
                self.log_group.name,
                self.log_stream_name,
                log_events,
            )
        elif service == "firehose":
            from moto.firehose import firehose_backends

            firehose_backends[self.account_id][self.region].send_log_event(
                destination_arn,
                filter_name,
                self.log_group.name,
                self.log_stream_name,
                log_events,
            )
        elif service == "kinesis":
            from moto.kinesis import kinesis_backends

            kinesis = kinesis_backends[self.account_id][self.region]
            kinesis.send_log_event(
                destination_arn,
                filter_name,
                self.log_group.name,
                self.log_stream_name,
                log_events,
            )

    def get_log_events(
        self,
        start_time: str,
        end_time: str,
        limit: int,
        next_token: Optional[str],
        start_from_head: str,
    ) -> Tuple[List[Dict[str, Any]], Optional[str], Optional[str]]:
        if limit is None:
            limit = 10000

        def filter_func(event: LogEvent) -> bool:
            if start_time and event.timestamp < start_time:
                return False

            if end_time and event.timestamp > end_time:
                return False

            return True

        def get_index_and_direction_from_token(
            token: Optional[str],
        ) -> Tuple[Optional[str], int]:
            if token is not None:
                try:
                    return token[0], int(token[2:])
                except Exception:
                    raise InvalidParameterException(
                        "The specified nextToken is invalid."
                    )
            return None, 0

        events = sorted(
            filter(filter_func, self.events), key=lambda event: event.timestamp
        )

        direction, index = get_index_and_direction_from_token(next_token)
        limit_index = limit - 1
        final_index = len(events) - 1

        if direction is None:
            if start_from_head:
                start_index = 0
                end_index = start_index + limit_index
            else:
                end_index = final_index
                start_index = end_index - limit_index
        elif direction == "f":
            start_index = index + 1
            end_index = start_index + limit_index
        elif direction == "b":
            end_index = index - 1
            start_index = end_index - limit_index
        else:
            raise InvalidParameterException("The specified nextToken is invalid.")

        if start_index < 0:
            start_index = 0
        elif start_index > final_index:
            return ([], f"b/{final_index:056d}", f"f/{final_index:056d}")

        if end_index > final_index:
            end_index = final_index
        elif end_index < 0:
            return ([], f"b/{0:056d}", f"f/{0:056d}")

        events_page = [
            event.to_response_dict() for event in events[start_index : end_index + 1]
        ]

        return (events_page, f"b/{start_index:056d}", f"f/{end_index:056d}")

    def filter_log_events(
        self, start_time: int, end_time: int, filter_pattern: str
    ) -> List[Dict[str, Any]]:
        def filter_func(event: LogEvent) -> bool:
            if start_time and event.timestamp < start_time:
                return False

            if end_time and event.timestamp > end_time:
                return False

            if not EventMessageFilter(filter_pattern).matches(event.message):
                return False

            return True

        events: List[Dict[str, Any]] = []
        for event in sorted(
            filter(filter_func, self.events), key=lambda x: x.timestamp
        ):
            event_obj = event.to_filter_dict()
            event_obj["logStreamName"] = self.log_stream_name
            events.append(event_obj)
        return events


class SubscriptionFilter(BaseModel):
    def __init__(
        self,
        name: str,
        log_group_name: str,
        filter_pattern: str,
        destination_arn: str,
        role_arn: str,
    ):
        self.name = name
        self.log_group_name = log_group_name
        self.filter_pattern = filter_pattern
        self.destination_arn = destination_arn
        self.role_arn = role_arn
        self.creation_time = int(unix_time_millis())

    def update(self, filter_pattern: str, destination_arn: str, role_arn: str) -> None:
        self.filter_pattern = filter_pattern
        self.destination_arn = destination_arn
        self.role_arn = role_arn

    def to_json(self) -> Dict[str, Any]:
        return {
            "filterName": self.name,
            "logGroupName": self.log_group_name,
            "filterPattern": self.filter_pattern,
            "destinationArn": self.destination_arn,
            "roleArn": self.role_arn,
            "distribution": "ByLogStream",
            "creationTime": self.creation_time,
        }


class LogGroup(CloudFormationModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        name: str,
        **kwargs: Any,
    ):
        self.name = name
        self.account_id = account_id
        self.region = region
        self.arn = (
            f"arn:{get_partition(region)}:logs:{region}:{account_id}:log-group:{name}"
        )
        self.creation_time = int(unix_time_millis())
        self.streams: Dict[str, LogStream] = dict()  # {name: LogStream}
        # AWS defaults to Never Expire for log group retention
        self.retention_in_days = kwargs.get("RetentionInDays")
        self.subscription_filters: Dict[str, SubscriptionFilter] = {}

        # The Amazon Resource Name (ARN) of the CMK to use when encrypting log data. It is optional.
        # Docs:
        # https://docs.aws.amazon.com/AmazonCloudWatchLogs/latest/APIReference/API_CreateLogGroup.html
        self.kms_key_id = kwargs.get("kmsKeyId")

    @staticmethod
    def cloudformation_name_type() -> str:
        return "LogGroupName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-logs-loggroup.html
        return "AWS::Logs::LogGroup"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "LogGroup":
        properties = cloudformation_json["Properties"]
        tags = properties.get("Tags", [])
        tags = dict([tag.values() for tag in tags])

        return logs_backends[account_id][region_name].create_log_group(
            resource_name, tags, **properties
        )

    def delete(self, account_id: str, region_name: str) -> None:
        backend = logs_backends[account_id][region_name]
        backend.delete_log_group(self.name)

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return self.arn
        raise UnformattedGetAttTemplateException()

    @property
    def physical_resource_id(self) -> str:
        return self.name

    def create_log_stream(self, log_stream_name: str) -> None:
        if log_stream_name in self.streams:
            raise ResourceAlreadyExistsException()
        stream = LogStream(log_group=self, name=log_stream_name)

        self.streams[log_stream_name] = stream

    def delete_log_stream(self, log_stream_name: str) -> None:
        if log_stream_name not in self.streams:
            raise ResourceNotFoundException()
        del self.streams[log_stream_name]

    def describe_log_streams(
        self,
        descending: bool,
        log_group_name: str,
        log_stream_name_prefix: str,
        order_by: str,
        limit: int,
        next_token: Optional[str] = None,
    ) -> Tuple[List[Dict[str, Any]], Optional[str]]:
        # responses only log_stream_name, creation_time, arn, stored_bytes when no events are stored.

        log_streams = [
            (name, stream.to_describe_dict())
            for name, stream in self.streams.items()
            if name.startswith(log_stream_name_prefix)
        ]

        def sorter(item: Any) -> Any:
            return (
                item[0]
                if order_by == "LogStreamName"
                else item[1].get("lastEventTimestamp", 0)
            )

        log_streams = sorted(log_streams, key=sorter, reverse=descending)
        first_index = 0
        if next_token:
            try:
                group, stream = next_token.split("@")
                if group != log_group_name:
                    raise ValueError()
                first_index = (
                    next(
                        index
                        for (index, e) in enumerate(log_streams)
                        if e[1]["logStreamName"] == stream
                    )
                    + 1
                )
            except (ValueError, StopIteration):
                first_index = 0
                log_streams = []

        last_index = first_index + limit
        if last_index > len(log_streams):
            last_index = len(log_streams)
        log_streams_page = [x[1] for x in log_streams[first_index:last_index]]
        new_token = None
        if log_streams_page and last_index < len(log_streams):
            new_token = f"{log_group_name}@{log_streams_page[-1]['logStreamName']}"

        return log_streams_page, new_token

    def put_log_events(
        self,
        log_stream_name: str,
        log_events: List[Dict[str, Any]],
    ) -> str:
        if log_stream_name not in self.streams:
            raise ResourceNotFoundException("The specified log stream does not exist.")
        stream = self.streams[log_stream_name]
        return stream.put_log_events(log_events)

    def get_log_events(
        self,
        log_stream_name: str,
        start_time: str,
        end_time: str,
        limit: int,
        next_token: Optional[str],
        start_from_head: str,
    ) -> Tuple[List[Dict[str, Any]], Optional[str], Optional[str]]:
        if log_stream_name not in self.streams:
            raise ResourceNotFoundException()
        stream = self.streams[log_stream_name]
        return stream.get_log_events(
            start_time,
            end_time,
            limit,
            next_token,
            start_from_head,
        )

    def filter_log_events(
        self,
        log_group_name: str,
        log_stream_names: List[str],
        start_time: int,
        end_time: int,
        limit: Optional[int],
        next_token: Optional[str],
        filter_pattern: str,
        interleaved: bool,
    ) -> Tuple[List[Dict[str, Any]], Optional[str], List[Dict[str, Any]]]:
        if not limit:
            limit = 10000
        streams = [
            stream
            for name, stream in self.streams.items()
            if not log_stream_names or name in log_stream_names
        ]

        events = []
        for stream in streams:
            events += stream.filter_log_events(start_time, end_time, filter_pattern)

        if interleaved:
            events = sorted(events, key=lambda event: event["timestamp"])

        first_index = 0
        if next_token:
            try:
                group, stream_name, event_id = next_token.split("@")
                if group != log_group_name:
                    raise ValueError()
                first_index = (
                    next(
                        index
                        for (index, e) in enumerate(events)
                        if e["logStreamName"] == stream_name
                        and e["eventId"] == event_id
                    )
                    + 1
                )
            except (ValueError, StopIteration):
                first_index = 0
                # AWS returns an empty list if it receives an invalid token.
                events = []

        last_index = first_index + limit
        if last_index > len(events):
            last_index = len(events)
        events_page = events[first_index:last_index]
        next_token = None
        if events_page and last_index < len(events):
            last_event = events_page[-1]
            next_token = f"{log_group_name}@{last_event['logStreamName']}@{last_event['eventId']}"

        searched_streams = [
            {"logStreamName": stream.log_stream_name, "searchedCompletely": True}
            for stream in streams
        ]
        return events_page, next_token, searched_streams

    def to_describe_dict(self) -> Dict[str, Any]:
        log_group = {
            "arn": self.arn,
            "creationTime": self.creation_time,
            "logGroupName": self.name,
            "metricFilterCount": 0,
            "storedBytes": sum(s.stored_bytes for s in self.streams.values()),
        }
        # AWS only returns retentionInDays if a value is set for the log group (ie. not Never Expire)
        if self.retention_in_days:
            log_group["retentionInDays"] = self.retention_in_days
        if self.kms_key_id:
            log_group["kmsKeyId"] = self.kms_key_id
        return log_group

    def set_retention_policy(self, retention_in_days: Optional[str]) -> None:
        self.retention_in_days = retention_in_days

    def describe_subscription_filters(self) -> Iterable[SubscriptionFilter]:
        return self.subscription_filters.values()

    def put_subscription_filter(
        self, filter_name: str, filter_pattern: str, destination_arn: str, role_arn: str
    ) -> None:
        # only two subscription filters can be associated with a log group
        if len(self.subscription_filters) == 2:
            raise LimitExceededException()

        # Update existing filter
        if filter_name in self.subscription_filters:
            self.subscription_filters[filter_name].update(
                filter_pattern, destination_arn, role_arn
            )
            return

        self.subscription_filters[filter_name] = SubscriptionFilter(
            name=filter_name,
            log_group_name=self.name,
            filter_pattern=filter_pattern,
            destination_arn=destination_arn,
            role_arn=role_arn,
        )

    def delete_subscription_filter(self, filter_name: str) -> None:
        if filter_name not in self.subscription_filters:
            raise ResourceNotFoundException(
                "The specified subscription filter does not exist."
            )

        self.subscription_filters.pop(filter_name)


class LogResourcePolicy(CloudFormationModel):
    def __init__(self, policy_name: str, policy_document: str):
        self.policy_name = policy_name
        self.policy_document = policy_document
        self.last_updated_time = int(unix_time_millis())

    def update(self, policy_document: str) -> None:
        self.policy_document = policy_document
        self.last_updated_time = int(unix_time_millis())

    def describe(self) -> Dict[str, Any]:
        return {
            "policyName": self.policy_name,
            "policyDocument": self.policy_document,
            "lastUpdatedTime": self.last_updated_time,
        }

    @property
    def physical_resource_id(self) -> str:
        return self.policy_name

    @staticmethod
    def cloudformation_name_type() -> str:
        return "PolicyName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-logs-resourcepolicy.html
        return "AWS::Logs::ResourcePolicy"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "LogResourcePolicy":
        properties = cloudformation_json["Properties"]
        policy_name = properties["PolicyName"]
        policy_document = properties["PolicyDocument"]
        return logs_backends[account_id][region_name].put_resource_policy(
            policy_name, policy_document
        )

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "LogResourcePolicy":
        properties = cloudformation_json["Properties"]
        policy_name = properties["PolicyName"]
        policy_document = properties["PolicyDocument"]

        backend = logs_backends[account_id][region_name]
        updated = backend.put_resource_policy(policy_name, policy_document)
        # TODO: move `update by replacement logic` to cloudformation. this is required for implementing rollbacks
        if original_resource.policy_name != policy_name:
            backend.delete_resource_policy(original_resource.policy_name)
        return updated

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        logs_backends[account_id][region_name].delete_resource_policy(resource_name)


class ExportTask(BaseModel):
    def __init__(
        self,
        task_id: str,
        task_name: str,
        log_group_name: str,
        destination: str,
        destination_prefix: str,
        from_time: int,
        to: int,
    ):
        self.task_id = task_id
        self.task_name = task_name
        self.log_group_name = log_group_name
        self.destination = destination
        self.destination_prefix = destination_prefix
        self.from_time = from_time
        self.to = to
        self.status = {"code": "active", "message": "Task is active"}

    def to_json(self) -> Dict[str, Any]:
        return {
            "taskId": self.task_id,
            "taskName": self.task_name,
            "logGroupName": self.log_group_name,
            "destination": self.destination,
            "destinationPrefix": self.destination_prefix,
            "from": self.from_time,
            "to": self.to,
            "status": self.status,
        }


class LogsBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.groups: Dict[str, LogGroup] = dict()
        self.filters = MetricFilters()
        self.queries: Dict[str, LogQuery] = dict()
        self.resource_policies: Dict[str, LogResourcePolicy] = dict()
        self.destinations: Dict[str, Destination] = dict()
        self.tagger = TaggingService()
        self.export_tasks: Dict[str, ExportTask] = dict()

    def create_log_group(
        self, log_group_name: str, tags: Dict[str, str], **kwargs: Any
    ) -> LogGroup:
        if log_group_name in self.groups:
            raise ResourceAlreadyExistsException()
        if len(log_group_name) > 512:
            raise InvalidParameterException(
                constraint="Member must have length less than or equal to 512",
                parameter="logGroupName",
                value=log_group_name,
            )
        self.groups[log_group_name] = LogGroup(
            self.account_id, self.region_name, log_group_name, **kwargs
        )
        self.tag_resource(self.groups[log_group_name].arn, tags)
        return self.groups[log_group_name]

    def ensure_log_group(self, log_group_name: str) -> None:
        if log_group_name in self.groups:
            return
        self.groups[log_group_name] = LogGroup(
            self.account_id,
            self.region_name,
            log_group_name,
        )

    def delete_log_group(self, log_group_name: str) -> None:
        if log_group_name not in self.groups:
            raise ResourceNotFoundException()
        del self.groups[log_group_name]

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore[misc]
    def describe_log_groups(
        self, log_group_name_prefix: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        groups = [
            group.to_describe_dict()
            for name, group in self.groups.items()
            if name.startswith(log_group_name_prefix or "")
        ]
        groups = sorted(groups, key=lambda x: x["logGroupName"])

        return groups

    def get_destination(self, destination_name: str) -> Destination:
        for destination in self.destinations:
            if self.destinations[destination].destination_name == destination_name:
                return self.destinations[destination]
        raise ResourceNotFoundException()

    def put_destination(
        self,
        destination_name: str,
        role_arn: str,
        target_arn: str,
        tags: Dict[str, str],
    ) -> Destination:
        for _, destination in self.destinations.items():
            if destination.destination_name == destination_name:
                if role_arn:
                    destination.role_arn = role_arn
                if target_arn:
                    destination.target_arn = target_arn
                return destination
        destination = Destination(
            self.account_id, self.region_name, destination_name, role_arn, target_arn
        )
        self.destinations[destination.arn] = destination
        self.tag_resource(destination.arn, tags)
        return destination

    def delete_destination(self, destination_name: str) -> None:
        destination = self.get_destination(destination_name)
        self.destinations.pop(destination.arn)
        return

    def describe_destinations(
        self, destination_name_prefix: str, limit: int, next_token: Optional[int] = None
    ) -> Tuple[List[Dict[str, Any]], Optional[int]]:
        if limit > 50:
            raise InvalidParameterException(
                constraint="Member must have value less than or equal to 50",
                parameter="limit",
                value=limit,
            )

        result = []
        for destination in self.destinations:
            result.append(self.destinations[destination].to_dict())
        if next_token:
            result = result[: int(next_token)]
        result = [
            destination
            for destination in result
            if destination["destinationName"].startswith(destination_name_prefix)
        ]
        return result, next_token

    def put_destination_policy(self, destination_name: str, access_policy: str) -> None:
        destination = self.get_destination(destination_name)
        destination.access_policy = access_policy
        return

    def create_log_stream(self, log_group_name: str, log_stream_name: str) -> None:
        if log_group_name not in self.groups:
            raise ResourceNotFoundException()
        log_group = self.groups[log_group_name]
        log_group.create_log_stream(log_stream_name)

    def ensure_log_stream(self, log_group_name: str, log_stream_name: str) -> None:
        if log_group_name not in self.groups:
            raise ResourceNotFoundException()

        if log_stream_name in self.groups[log_group_name].streams:
            return

        self.create_log_stream(log_group_name, log_stream_name)

    def delete_log_stream(self, log_group_name: str, log_stream_name: str) -> None:
        if log_group_name not in self.groups:
            raise ResourceNotFoundException()
        log_group = self.groups[log_group_name]
        log_group.delete_log_stream(log_stream_name)

    def describe_log_streams(
        self,
        descending: bool,
        limit: int,
        log_group_name: str,
        log_stream_name_prefix: str,
        next_token: Optional[str],
        order_by: str,
    ) -> Tuple[List[Dict[str, Any]], Optional[str]]:
        if log_group_name not in self.groups:
            raise ResourceNotFoundException()
        if limit > 50:
            raise InvalidParameterException(
                constraint="Member must have value less than or equal to 50",
                parameter="limit",
                value=limit,
            )
        if order_by not in ["LogStreamName", "LastEventTime"]:
            raise InvalidParameterException(
                constraint="Member must satisfy enum value set: [LogStreamName, LastEventTime]",
                parameter="orderBy",
                value=order_by,
            )
        if order_by == "LastEventTime" and log_stream_name_prefix:
            raise InvalidParameterException(
                msg="Cannot order by LastEventTime with a logStreamNamePrefix."
            )
        log_group = self.groups[log_group_name]
        return log_group.describe_log_streams(
            descending=descending,
            limit=limit,
            log_group_name=log_group_name,
            log_stream_name_prefix=log_stream_name_prefix,
            next_token=next_token,
            order_by=order_by,
        )

    def put_log_events(
        self,
        log_group_name: str,
        log_stream_name: str,
        log_events: List[Dict[str, Any]],
    ) -> Tuple[str, Dict[str, Any]]:
        """
        The SequenceToken-parameter is not yet implemented
        """
        if log_group_name not in self.groups:
            raise ResourceNotFoundException()
        log_group = self.groups[log_group_name]

        # Only events from the last 14 days or 2 hours in the future are accepted
        rejected_info = {}
        allowed_events = []
        last_timestamp = None
        oldest = int(unix_time_millis(utcnow() - timedelta(days=14)))
        newest = int(unix_time_millis(utcnow() + timedelta(hours=2)))
        for idx, event in enumerate(log_events):
            if last_timestamp and last_timestamp > event["timestamp"]:
                raise InvalidParameterException(
                    "Log events in a single PutLogEvents request must be in chronological order."
                )
            if event["timestamp"] < oldest:
                rejected_info["tooOldLogEventEndIndex"] = idx
            elif event["timestamp"] > newest:
                rejected_info["tooNewLogEventStartIndex"] = idx
            else:
                allowed_events.append(event)
            last_timestamp = event["timestamp"]

        token = log_group.put_log_events(log_stream_name, allowed_events)
        return token, rejected_info

    def get_log_events(
        self,
        log_group_name: str,
        log_stream_name: str,
        start_time: str,
        end_time: str,
        limit: int,
        next_token: Optional[str],
        start_from_head: str,
    ) -> Tuple[List[Dict[str, Any]], Optional[str], Optional[str]]:
        if log_group_name not in self.groups:
            raise ResourceNotFoundException()
        if limit and limit > 10000:
            raise InvalidParameterException(
                constraint="Member must have value less than or equal to 10000",
                parameter="limit",
                value=limit,
            )
        log_group = self.groups[log_group_name]
        return log_group.get_log_events(
            log_stream_name, start_time, end_time, limit, next_token, start_from_head
        )

    def filter_log_events(
        self,
        log_group_name: str,
        log_stream_names: List[str],
        start_time: int,
        end_time: int,
        limit: Optional[int],
        next_token: Optional[str],
        filter_pattern: str,
        interleaved: bool,
    ) -> Tuple[List[Dict[str, Any]], Optional[str], List[Dict[str, Any]]]:
        """
        The following filter patterns are currently supported: Single Terms, Multiple Terms, Exact Phrases.
        If the pattern is not supported, all events are returned.
        """
        if log_group_name not in self.groups:
            raise ResourceNotFoundException()
        if limit and limit > 10000:
            raise InvalidParameterException(
                constraint="Member must have value less than or equal to 10000",
                parameter="limit",
                value=limit,
            )
        log_group = self.groups[log_group_name]
        return log_group.filter_log_events(
            log_group_name,
            log_stream_names,
            start_time,
            end_time,
            limit,
            next_token,
            filter_pattern,
            interleaved,
        )

    def put_retention_policy(self, log_group_name: str, retention_in_days: str) -> None:
        if log_group_name not in self.groups:
            raise ResourceNotFoundException()
        self.groups[log_group_name].set_retention_policy(retention_in_days)

    def delete_retention_policy(self, log_group_name: str) -> None:
        if log_group_name not in self.groups:
            raise ResourceNotFoundException()
        self.groups[log_group_name].set_retention_policy(None)

    def describe_resource_policies(self) -> List[LogResourcePolicy]:
        """
        Return list of resource policies.

        The next_token and limit arguments are ignored.  The maximum
        number of resource policies per region is a small number (less
        than 50), so pagination isn't needed.
        """

        return list(self.resource_policies.values())

    def put_resource_policy(
        self, policy_name: str, policy_doc: str
    ) -> LogResourcePolicy:
        """
        Creates/updates resource policy and return policy object
        """
        if policy_name in self.resource_policies:
            policy = self.resource_policies[policy_name]
            policy.update(policy_doc)
            return policy
        if len(self.resource_policies) == MAX_RESOURCE_POLICIES_PER_REGION:
            raise LimitExceededException()
        policy = LogResourcePolicy(policy_name, policy_doc)
        self.resource_policies[policy_name] = policy
        return policy

    def delete_resource_policy(self, policy_name: str) -> None:
        """
        Remove resource policy with a policy name matching given name.
        """
        if policy_name not in self.resource_policies:
            raise ResourceNotFoundException(
                msg=f"Policy with name [{policy_name}] does not exist"
            )
        del self.resource_policies[policy_name]

    def list_tags_log_group(self, log_group_name: str) -> Dict[str, str]:
        if log_group_name not in self.groups:
            raise ResourceNotFoundException()
        log_group = self.groups[log_group_name]
        return self.list_tags_for_resource(log_group.arn)

    def tag_log_group(self, log_group_name: str, tags: Dict[str, str]) -> None:
        if log_group_name not in self.groups:
            raise ResourceNotFoundException()
        log_group = self.groups[log_group_name]
        self.tag_resource(log_group.arn, tags)

    def untag_log_group(self, log_group_name: str, tags: List[str]) -> None:
        if log_group_name not in self.groups:
            raise ResourceNotFoundException()
        log_group = self.groups[log_group_name]
        self.untag_resource(log_group.arn, tags)

    def put_metric_filter(
        self,
        filter_name: str,
        filter_pattern: str,
        log_group_name: str,
        metric_transformations: str,
    ) -> None:
        self.filters.add_filter(
            filter_name, filter_pattern, log_group_name, metric_transformations
        )

    def describe_metric_filters(
        self,
        prefix: Optional[str] = None,
        log_group_name: Optional[str] = None,
        metric_name: Optional[str] = None,
        metric_namespace: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        filters = self.filters.get_matching_filters(
            prefix, log_group_name, metric_name, metric_namespace
        )
        return filters

    def delete_metric_filter(
        self, filter_name: Optional[str] = None, log_group_name: Optional[str] = None
    ) -> None:
        self.filters.delete_filter(filter_name, log_group_name)

    def describe_subscription_filters(
        self, log_group_name: str
    ) -> Iterable[SubscriptionFilter]:
        log_group = self.groups.get(log_group_name)

        if not log_group:
            raise ResourceNotFoundException()

        return log_group.describe_subscription_filters()

    def put_subscription_filter(
        self,
        log_group_name: str,
        filter_name: str,
        filter_pattern: str,
        destination_arn: str,
        role_arn: str,
    ) -> None:
        log_group = self.groups.get(log_group_name)

        if not log_group:
            raise ResourceNotFoundException()

        service = destination_arn.split(":")[2]
        if service == "lambda":
            from moto.awslambda.utils import get_backend

            try:
                get_backend(self.account_id, self.region_name).get_function(
                    destination_arn
                )
            # no specific permission check implemented
            except Exception:
                raise InvalidParameterException(
                    "Could not execute the lambda function. Make sure you "
                    "have given CloudWatch Logs permission to execute your "
                    "function."
                )
        elif service == "firehose":
            from moto.firehose import firehose_backends

            firehose = firehose_backends[self.account_id][
                self.region_name
            ].lookup_name_from_arn(destination_arn)
            if not firehose:
                raise InvalidParameterException(
                    "Could not deliver test message to specified Firehose "
                    "stream. Check if the given Firehose stream is in ACTIVE "
                    "state."
                )
        elif service == "kinesis":
            from moto.kinesis import kinesis_backends

            kinesis = kinesis_backends[self.account_id][self.region_name]
            try:
                kinesis.describe_stream(stream_arn=destination_arn, stream_name=None)
            except Exception:
                raise InvalidParameterException(
                    "Could not deliver test message to specified Kinesis stream. Verify the stream exists "
                )
        else:
            # TODO: support Kinesis stream destinations
            raise InvalidParameterException(
                f"Service '{service}' has not implemented for "
                f"put_subscription_filter()"
            )

        log_group.put_subscription_filter(
            filter_name, filter_pattern, destination_arn, role_arn
        )

    def delete_subscription_filter(self, log_group_name: str, filter_name: str) -> None:
        log_group = self.groups.get(log_group_name)

        if not log_group:
            raise ResourceNotFoundException()

        log_group.delete_subscription_filter(filter_name)

    def start_query(
        self,
        log_group_names: List[str],
        start_time: int,
        end_time: int,
        query_string: str,
    ) -> str:
        for log_group_name in log_group_names:
            if log_group_name not in self.groups:
                raise ResourceNotFoundException()
        log_groups = [self.groups[name] for name in log_group_names]

        query_id = str(mock_random.uuid1())
        self.queries[query_id] = LogQuery(
            query_id, start_time, end_time, query_string, log_groups
        )
        return query_id

    def describe_queries(
        self, log_stream_name: str, status: Optional[str]
    ) -> List[LogQuery]:
        """
        Pagination is not yet implemented
        """
        queries: List[LogQuery] = []
        for query in self.queries.values():
            if log_stream_name in query.log_group_names and (
                not status or status == query.status
            ):
                queries.append(query)
        return queries

    def get_query_results(self, query_id: str) -> LogQuery:
        """
        Not all query commands are implemented yet. Please raise an issue if you encounter unexpected results.
        """
        return self.queries[query_id]

    def create_export_task(
        self,
        taskName: str,
        logGroupName: str,
        destination: str,
        destinationPrefix: str,
        fromTime: int,
        to: int,
    ) -> str:
        try:
            s3_backends[self.account_id][self.partition].get_bucket(destination)
        except MissingBucket:
            raise InvalidParameterException(
                "The given bucket does not exist. Please make sure the bucket is valid."
            )
        if logGroupName not in self.groups:
            raise ResourceNotFoundException()
        task_id = str(mock_random.uuid4())
        self.export_tasks[task_id] = ExportTask(
            task_id,
            taskName,
            logGroupName,
            destination,
            destinationPrefix,
            fromTime,
            to,
        )

        s3_backends[self.account_id][self.partition].put_object(
            bucket_name=destination,
            key_name="aws-logs-write-test",
            value=b"Permission Check Successful",
        )

        if fromTime <= to:
            for stream_name in self.groups[logGroupName].streams.keys():
                logs, _, _ = self.filter_log_events(
                    log_group_name=logGroupName,
                    log_stream_names=[stream_name],
                    start_time=fromTime,
                    end_time=to,
                    limit=None,
                    next_token=None,
                    filter_pattern="",
                    interleaved=False,
                )
                raw_logs = "\n".join(
                    [
                        f"{datetime.fromtimestamp(log['timestamp']/1000).strftime('%Y-%m-%dT%H:%M:%S.000Z')} {log['message']}"
                        for log in logs
                    ]
                )
                folder = str(mock_random.uuid4()) + "/" + stream_name.replace("/", "-")
                key_name = f"{destinationPrefix}/{folder}/000000.gz"
                s3_backends[self.account_id][self.partition].put_object(
                    bucket_name=destination,
                    key_name=key_name,
                    value=gzip_compress(raw_logs.encode("utf-8")),
                )
            self.export_tasks[task_id].status["code"] = "COMPLETED"
            self.export_tasks[task_id].status["message"] = "Completed successfully"

        return task_id

    def describe_export_tasks(self, task_id: str) -> List[ExportTask]:
        """
        Pagination is not yet implemented
        """
        if task_id:
            if task_id not in self.export_tasks:
                raise ResourceNotFoundException()
            return [self.export_tasks[task_id]]
        else:
            return list(self.export_tasks.values())

    def list_tags_for_resource(self, resource_arn: str) -> Dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_arn)

    def tag_resource(self, arn: str, tags: Dict[str, str]) -> None:
        self.tagger.tag_resource(arn, TaggingService.convert_dict_to_tags_input(tags))

    def untag_resource(self, arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(arn, tag_keys)


logs_backends = BackendDict(LogsBackend, "logs")
