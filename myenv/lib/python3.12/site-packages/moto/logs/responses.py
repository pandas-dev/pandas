import json
import re
from typing import Any, Callable, Optional

from moto.core.responses import BaseResponse

from .exceptions import InvalidParameterException
from .models import LogsBackend, logs_backends

# See http://docs.aws.amazon.com/AmazonCloudWatchLogs/latest/APIReference/Welcome.html


REGEX_LOG_GROUP_NAME = r"[-._\/#A-Za-z0-9]+"


def validate_param(
    param_name: str,
    param_value: str,
    constraint: str,
    constraint_expression: Callable[[str], bool],
    pattern: Optional[str] = None,
) -> None:
    try:
        assert constraint_expression(param_value)
    except (AssertionError, TypeError):
        raise InvalidParameterException(
            constraint=constraint, parameter=param_name, value=param_value
        )
    if pattern and param_value:
        try:
            assert re.fullmatch(pattern, param_value)
        except (AssertionError, TypeError):
            raise InvalidParameterException(
                constraint=f"Must match pattern: {pattern}",
                parameter=param_name,
                value=param_value,
            )


class LogsResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="logs")

    @property
    def logs_backend(self) -> LogsBackend:
        return logs_backends[self.current_account][self.region]

    def _get_validated_param(
        self,
        param: str,
        constraint: str,
        constraint_expression: Callable[[str], bool],
        pattern: Optional[str] = None,
    ) -> Any:
        param_value = self._get_param(param)
        validate_param(param, param_value, constraint, constraint_expression, pattern)
        return param_value

    def put_metric_filter(self) -> str:
        filter_name = self._get_validated_param(
            "filterName",
            "Minimum length of 1. Maximum length of 512.",
            lambda x: 1 <= len(x) <= 512,
            pattern="[^:*]*",
        )
        filter_pattern = self._get_validated_param(
            "filterPattern",
            "Minimum length of 0. Maximum length of 1024.",
            lambda x: 0 <= len(x) <= 1024,
        )
        log_group_name = self._get_validated_param(
            "logGroupName",
            "Minimum length of 1. Maximum length of 512.",
            lambda x: 1 <= len(x) <= 512,
            pattern=REGEX_LOG_GROUP_NAME,
        )
        metric_transformations = self._get_validated_param(
            "metricTransformations", "Fixed number of 1 item.", lambda x: len(x) == 1
        )

        self.logs_backend.put_metric_filter(
            filter_name, filter_pattern, log_group_name, metric_transformations
        )

        return ""

    def describe_metric_filters(self) -> str:
        filter_name_prefix = self._get_validated_param(
            "filterNamePrefix",
            "Minimum length of 1. Maximum length of 512.",
            lambda x: x is None or 1 <= len(x) <= 512,
            pattern="[^:*]*",
        )
        log_group_name = self._get_validated_param(
            "logGroupName",
            "Minimum length of 1. Maximum length of 512",
            lambda x: x is None or 1 <= len(x) <= 512,
            pattern=REGEX_LOG_GROUP_NAME,
        )
        metric_name = self._get_validated_param(
            "metricName",
            "Maximum length of 255.",
            lambda x: x is None or len(x) <= 255,
            pattern="[^:*$]*",
        )
        metric_namespace = self._get_validated_param(
            "metricNamespace",
            "Maximum length of 255.",
            lambda x: x is None or len(x) <= 255,
            pattern="[0-9A-Za-z\\.\\-_/#:]*",
        )
        next_token = self._get_validated_param(
            "nextToken", "Minimum length of 1.", lambda x: x is None or 1 <= len(x)
        )

        if metric_name and not metric_namespace:
            raise InvalidParameterException(
                constraint=f'{"If you include the metricName parameter in your request, "}'
                f'{"you must also include the metricNamespace parameter."}',
                parameter="metricNamespace",
                value=metric_namespace,
            )
        if metric_namespace and not metric_name:
            raise InvalidParameterException(
                constraint=f'{"If you include the metricNamespace parameter in your request, "}'
                f'{"you must also include the metricName parameter."}',
                parameter="metricName",
                value=metric_name,
            )

        filters = self.logs_backend.describe_metric_filters(
            filter_name_prefix, log_group_name, metric_name, metric_namespace
        )
        return json.dumps({"metricFilters": filters, "nextToken": next_token})

    def delete_metric_filter(self) -> str:
        filter_name = self._get_validated_param(
            "filterName",
            "Minimum length of 1. Maximum length of 512.",
            lambda x: 1 <= len(x) <= 512,
            pattern="[^:*]*$",
        )
        log_group_name = self._get_validated_param(
            "logGroupName",
            "Minimum length of 1. Maximum length of 512.",
            lambda x: 1 <= len(x) <= 512,
            pattern=REGEX_LOG_GROUP_NAME,
        )

        self.logs_backend.delete_metric_filter(filter_name, log_group_name)
        return ""

    def create_log_group(self) -> str:
        log_group_name = self._get_param("logGroupName")
        tags = self._get_param("tags")
        kms_key_id = self._get_param("kmsKeyId")

        self.logs_backend.create_log_group(log_group_name, tags, kmsKeyId=kms_key_id)
        return ""

    def delete_log_group(self) -> str:
        log_group_name = self._get_param("logGroupName")
        self.logs_backend.delete_log_group(log_group_name)
        return ""

    def describe_log_groups(self) -> str:
        log_group_name_prefix = self._get_param("logGroupNamePrefix")
        next_token = self._get_param("nextToken")
        limit = self._get_param("limit", 50)
        if limit > 50:
            raise InvalidParameterException(
                constraint="Member must have value less than or equal to 50",
                parameter="limit",
                value=limit,
            )
        groups, next_token = self.logs_backend.describe_log_groups(
            limit=limit,
            log_group_name_prefix=log_group_name_prefix,
            next_token=next_token,
        )
        result = {"logGroups": groups}
        if next_token:
            result["nextToken"] = next_token
        return json.dumps(result)

    def put_destination(self) -> str:
        destination_name = self._get_param("destinationName")
        role_arn = self._get_param("roleArn")
        target_arn = self._get_param("targetArn")
        tags = self._get_param("tags")

        destination = self.logs_backend.put_destination(
            destination_name,
            role_arn,
            target_arn,
            tags,
        )
        result = {"destination": destination.to_dict()}
        return json.dumps(result)

    def delete_destination(self) -> str:
        destination_name = self._get_param("destinationName")
        self.logs_backend.delete_destination(destination_name)
        return ""

    def describe_destinations(self) -> str:
        destination_name_prefix = self._get_param("DestinationNamePrefix")
        limit = self._get_param("limit", 50)
        next_token = self._get_param("nextToken")

        destinations, next_token = self.logs_backend.describe_destinations(
            destination_name_prefix, int(limit), next_token
        )

        result = {"destinations": destinations, "nextToken": next_token}
        return json.dumps(result)

    def put_destination_policy(self) -> str:
        access_policy = self._get_param("accessPolicy")
        destination_name = self._get_param("destinationName")

        self.logs_backend.put_destination_policy(destination_name, access_policy)
        return ""

    def create_log_stream(self) -> str:
        log_group_name = self._get_param("logGroupName")
        log_stream_name = self._get_param("logStreamName")
        self.logs_backend.create_log_stream(log_group_name, log_stream_name)
        return ""

    def delete_log_stream(self) -> str:
        log_group_name = self._get_param("logGroupName")
        log_stream_name = self._get_param("logStreamName")
        self.logs_backend.delete_log_stream(log_group_name, log_stream_name)
        return ""

    def describe_log_streams(self) -> str:
        log_group_name = self._get_param("logGroupName")
        log_stream_name_prefix = self._get_param("logStreamNamePrefix", "")
        descending = self._get_param("descending", False)
        limit = self._get_param("limit", 50)
        next_token = self._get_param("nextToken")
        order_by = self._get_param("orderBy", "LogStreamName")

        streams, next_token = self.logs_backend.describe_log_streams(
            descending,
            limit,
            log_group_name,
            log_stream_name_prefix,
            next_token,
            order_by,
        )
        return json.dumps({"logStreams": streams, "nextToken": next_token})

    def put_log_events(self) -> str:
        log_group_name = self._get_param("logGroupName")
        log_stream_name = self._get_param("logStreamName")
        log_events = self._get_param("logEvents")

        next_sequence_token, rejected_info = self.logs_backend.put_log_events(
            log_group_name, log_stream_name, log_events
        )
        if rejected_info:
            return json.dumps(
                {
                    "nextSequenceToken": next_sequence_token,
                    "rejectedLogEventsInfo": rejected_info,
                }
            )
        else:
            return json.dumps({"nextSequenceToken": next_sequence_token})

    def get_log_events(self) -> str:
        log_group_name = self._get_param("logGroupName")
        log_stream_name = self._get_param("logStreamName")
        start_time = self._get_param("startTime")
        end_time = self._get_param("endTime")
        limit = self._get_param("limit")
        next_token = self._get_param("nextToken")
        start_from_head = self._get_param("startFromHead", False)

        (
            events,
            next_backward_token,
            next_forward_token,
        ) = self.logs_backend.get_log_events(
            log_group_name,
            log_stream_name,
            start_time,
            end_time,
            limit,
            next_token,
            start_from_head,
        )
        return json.dumps(
            {
                "events": events,
                "nextBackwardToken": next_backward_token,
                "nextForwardToken": next_forward_token,
            }
        )

    def filter_log_events(self) -> str:
        log_group_name = self._get_param("logGroupName")
        log_stream_names = self._get_param("logStreamNames", [])
        start_time = self._get_param("startTime")
        # impl, see: http://docs.aws.amazon.com/AmazonCloudWatch/latest/logs/FilterAndPatternSyntax.html
        filter_pattern = self._get_param("filterPattern")
        interleaved = self._get_param("interleaved", False)
        end_time = self._get_param("endTime")
        limit = self._get_param("limit")
        next_token = self._get_param("nextToken")

        events, next_token, searched_streams = self.logs_backend.filter_log_events(
            log_group_name,
            log_stream_names,
            start_time,
            end_time,
            limit,
            next_token,
            filter_pattern,
            interleaved,
        )
        return json.dumps(
            {
                "events": events,
                "nextToken": next_token,
                "searchedLogStreams": searched_streams,
            }
        )

    def put_retention_policy(self) -> str:
        log_group_name = self._get_param("logGroupName")
        retention_in_days = self._get_param("retentionInDays")
        self.logs_backend.put_retention_policy(log_group_name, retention_in_days)
        return ""

    def delete_retention_policy(self) -> str:
        log_group_name = self._get_param("logGroupName")
        self.logs_backend.delete_retention_policy(log_group_name)
        return ""

    def describe_resource_policies(self) -> str:
        policies = self.logs_backend.describe_resource_policies()
        return json.dumps({"resourcePolicies": [p.describe() for p in policies]})

    def put_resource_policy(self) -> str:
        policy_name = self._get_param("policyName")
        policy_doc = self._get_param("policyDocument")
        policy = self.logs_backend.put_resource_policy(policy_name, policy_doc)
        return json.dumps({"resourcePolicy": policy.describe()})

    def delete_resource_policy(self) -> str:
        policy_name = self._get_param("policyName")
        self.logs_backend.delete_resource_policy(policy_name)
        return ""

    def list_tags_log_group(self) -> str:
        log_group_name = self._get_param("logGroupName")
        tags = self.logs_backend.list_tags_log_group(log_group_name)
        return json.dumps({"tags": tags})

    def tag_log_group(self) -> str:
        log_group_name = self._get_param("logGroupName")
        tags = self._get_param("tags")
        self.logs_backend.tag_log_group(log_group_name, tags)
        return ""

    def untag_log_group(self) -> str:
        log_group_name = self._get_param("logGroupName")
        tags = self._get_param("tags")
        self.logs_backend.untag_log_group(log_group_name, tags)
        return ""

    def describe_subscription_filters(self) -> str:
        log_group_name = self._get_param("logGroupName")

        _filters = self.logs_backend.describe_subscription_filters(log_group_name)

        return json.dumps({"subscriptionFilters": [f.to_json() for f in _filters]})

    def put_subscription_filter(self) -> str:
        log_group_name = self._get_param("logGroupName")
        filter_name = self._get_param("filterName")
        filter_pattern = self._get_param("filterPattern")
        destination_arn = self._get_param("destinationArn")
        role_arn = self._get_param("roleArn")

        self.logs_backend.put_subscription_filter(
            log_group_name, filter_name, filter_pattern, destination_arn, role_arn
        )

        return ""

    def delete_subscription_filter(self) -> str:
        log_group_name = self._get_param("logGroupName")
        filter_name = self._get_param("filterName")

        self.logs_backend.delete_subscription_filter(log_group_name, filter_name)

        return ""

    def start_query(self) -> str:
        log_group_name = self._get_param("logGroupName")
        log_group_names = self._get_param("logGroupNames")
        start_time = self._get_int_param("startTime")
        end_time = self._get_int_param("endTime")
        query_string = self._get_param("queryString")

        if log_group_name and log_group_names:
            raise InvalidParameterException()

        if log_group_name:
            log_group_names = [log_group_name]

        query_id = self.logs_backend.start_query(
            log_group_names, start_time, end_time, query_string
        )

        return json.dumps({"queryId": f"{query_id}"})

    def describe_queries(self) -> str:
        log_group_name = self._get_param("logGroupName")
        status = self._get_param("status")
        queries = self.logs_backend.describe_queries(log_group_name, status)
        return json.dumps(
            {"queries": [query.to_json(log_group_name) for query in queries]}
        )

    def get_query_results(self) -> str:
        query_id = self._get_param("queryId")
        query = self.logs_backend.get_query_results(query_id)
        return json.dumps(query.to_result_json())

    def create_export_task(self) -> str:
        task_id = self.logs_backend.create_export_task(
            logGroupName=self._get_param("logGroupName"),
            fromTime=self._get_int_param("from"),
            to=self._get_int_param("to"),
            destination=self._get_param("destination"),
            destinationPrefix=self._get_param("destinationPrefix", "exportedlogs"),
            taskName=self._get_param("taskName"),
        )
        return json.dumps(dict(taskId=str(task_id)))

    def describe_export_tasks(self) -> str:
        task_id = self._get_param("taskId")

        tasks = self.logs_backend.describe_export_tasks(task_id=task_id)
        return json.dumps({"exportTasks": [t.to_json() for t in tasks]})

    def list_tags_for_resource(self) -> str:
        resource_arn = self._get_param("resourceArn")
        tags = self.logs_backend.list_tags_for_resource(resource_arn)
        return json.dumps({"tags": tags})

    def tag_resource(self) -> str:
        resource_arn = self._get_param("resourceArn")
        tags = self._get_param("tags")
        self.logs_backend.tag_resource(resource_arn, tags)
        return "{}"

    def untag_resource(self) -> str:
        resource_arn = self._get_param("resourceArn")
        tag_keys = self._get_param("tagKeys")
        self.logs_backend.untag_resource(resource_arn, tag_keys)
        return "{}"
