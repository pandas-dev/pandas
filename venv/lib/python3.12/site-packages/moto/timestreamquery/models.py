"""TimestreamQueryBackend class with methods for supported APIs."""

from typing import Any, Dict, List, Optional, Union
from uuid import uuid4

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.utilities.utils import get_partition

from .exceptions import ResourceNotFound


class ScheduledQuery(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        query_string: str,
        schedule_configuration: Dict[str, str],
        notification_configuration: Dict[str, Dict[str, str]],
        target_configuration: Optional[Dict[str, Any]],
        scheduled_query_execution_role_arn: str,
        tags: Optional[List[Dict[str, str]]],
        kms_key_id: Optional[str],
        error_report_configuration: Optional[Dict[str, Dict[str, str]]],
    ):
        self.account_id = account_id
        self.region_name = region_name
        self.name = name
        self.query_string = query_string
        self.schedule_configuration = schedule_configuration
        self.notification_configuration = notification_configuration
        self.target_configuration = target_configuration
        self.scheduled_query_execution_role_arn = scheduled_query_execution_role_arn
        self.tags = tags
        self.kms_key_id = kms_key_id
        self.error_report_configuration = error_report_configuration

        self.created_on = unix_time()
        self.updated_on = unix_time()

        self.arn = f"arn:{get_partition(region_name)}:timestream:{region_name}:{account_id}:scheduled-query/{name}"
        self.state = "ENABLED"

    def description(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "Name": self.name,
            "CreationTime": self.created_on,
            "State": self.state,
            "QueryString": self.query_string,
            "ScheduleConfiguration": self.schedule_configuration,
            "NotificationConfiguration": self.notification_configuration,
            "TargetConfiguration": self.target_configuration,
            "ScheduledQueryExecutionRoleArn": self.scheduled_query_execution_role_arn,
            "KmsKeyId": self.kms_key_id,
            "ErrorReportConfiguration": self.error_report_configuration,
        }


class TimestreamQueryBackend(BaseBackend):
    """Implementation of TimestreamQuery APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.scheduled_queries: Dict[str, ScheduledQuery] = {}

        self.query_result_queue: Dict[Optional[str], List[Dict[str, Any]]] = {}
        self.query_results: Dict[str, Dict[str, Any]] = {}

    def create_scheduled_query(
        self,
        name: str,
        query_string: str,
        schedule_configuration: Dict[str, str],
        notification_configuration: Dict[str, Dict[str, str]],
        target_configuration: Optional[Dict[str, Any]],
        scheduled_query_execution_role_arn: str,
        tags: Optional[List[Dict[str, str]]],
        kms_key_id: Optional[str],
        error_report_configuration: Dict[str, Dict[str, str]],
    ) -> ScheduledQuery:
        query = ScheduledQuery(
            account_id=self.account_id,
            region_name=self.region_name,
            name=name,
            query_string=query_string,
            schedule_configuration=schedule_configuration,
            notification_configuration=notification_configuration,
            target_configuration=target_configuration,
            scheduled_query_execution_role_arn=scheduled_query_execution_role_arn,
            tags=tags,
            kms_key_id=kms_key_id,
            error_report_configuration=error_report_configuration,
        )
        self.scheduled_queries[query.arn] = query
        return query

    def delete_scheduled_query(self, scheduled_query_arn: str) -> None:
        self.scheduled_queries.pop(scheduled_query_arn, None)

    def describe_scheduled_query(self, scheduled_query_arn: str) -> ScheduledQuery:
        if scheduled_query_arn not in self.scheduled_queries:
            raise ResourceNotFound(scheduled_query_arn)
        return self.scheduled_queries[scheduled_query_arn]

    def update_scheduled_query(self, scheduled_query_arn: str, state: str) -> None:
        query = self.scheduled_queries[scheduled_query_arn]
        query.state = state

    def query(self, query_string: str) -> Dict[str, Any]:
        """
        Moto does not have a builtin time-series Database, so calling this endpoint will return zero results by default.

        You can use a dedicated API to configuring a queue of expected results.

        An example invocation looks like this:

        .. sourcecode:: python

            first_result = {
                'QueryId': 'some_id',
                'Rows': [...],
                'ColumnInfo': [...],
                'QueryStatus': ...
            }
            result_for_unknown_query_string = {
                'QueryId': 'unknown',
                'Rows': [...],
                'ColumnInfo': [...],
                'QueryStatus': ...
            }
            expected_results = {
                "account_id": "123456789012",  # This is the default - can be omitted
                "region": "us-east-1",  # This is the default - can be omitted
                "results": {
                    # Use the exact querystring, and a list of results for it
                    # For example
                    "SELECT data FROM mytable": [first_result, ...],
                    # Use None if the exact querystring is unknown/irrelevant
                    None: [result_for_unknown_query_string, ...],
                }
            }
            requests.post(
                "http://motoapi.amazonaws.com/moto-api/static/timestream/query-results",
                json=expected_results,
            )

        When calling `query(QueryString='SELECT data FROM mytable')`, the `first_result` will be returned.
        Call the query again for the second result, and so on.

        If you don't know the exact query strings, use the `None`-key. In the above example, when calling `SELECT something FROM unknown`, there are no results for that specific query, so `result_for_unknown_query_string` will be returned.

        Results for unknown queries are cached, so calling `SELECT something FROM unknown` will return the same result.

        """
        if self.query_result_queue.get(query_string):
            return self.query_result_queue[query_string].pop(0)
        if result := self.query_results.get(query_string):
            return result
        if self.query_result_queue.get(None):
            self.query_results[query_string] = self.query_result_queue[None].pop(0)
            return self.query_results[query_string]
        return {"QueryId": str(uuid4()), "Rows": [], "ColumnInfo": []}

    def describe_endpoints(self) -> List[Dict[str, Union[str, int]]]:
        # https://docs.aws.amazon.com/timestream/latest/developerguide/Using-API.endpoint-discovery.how-it-works.html
        # Usually, the address look like this:
        # query-cell1.timestream.us-east-1.amazonaws.com
        # Where 'cell1' can be any number, 'cell2', 'cell3', etc - whichever endpoint happens to be available for that particular account
        # We don't implement a cellular architecture in Moto though, so let's keep it simple
        return [
            {
                "Address": f"query.timestream.{self.region_name}.amazonaws.com",
                "CachePeriodInMinutes": 1440,
            }
        ]


timestreamquery_backends = BackendDict(
    TimestreamQueryBackend,
    "timestream-query",
    additional_regions=[
        "us-east-1",
        "us-east-2",
        "us-west-2",
        "eu-central-1",
        "eu-west-1",
        "ap-southeast-2",
        "ap-northeast-1",
    ],
)
