import re
from datetime import datetime
from typing import Any, Dict, Iterable, Iterator, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.utils import iso_8601_datetime_without_milliseconds
from moto.moto_api._internal import mock_random as random
from moto.redshiftdata.exceptions import ResourceNotFoundException, ValidationException


class Statement(Iterable[Tuple[str, Any]]):
    def __init__(
        self,
        cluster_identifier: str,
        database: str,
        db_user: str,
        query_parameters: List[Dict[str, str]],
        query_string: str,
        secret_arn: str,
    ):
        now = iso_8601_datetime_without_milliseconds(datetime.now())

        self.id = str(random.uuid4())
        self.cluster_identifier = cluster_identifier
        self.created_at = now
        self.database = database
        self.db_user = db_user
        self.duration = 0
        self.has_result_set = False
        self.query_parameters = query_parameters
        self.query_string = query_string
        self.redshift_pid = random.randint(0, 99999)
        self.redshift_query_id = random.randint(0, 99999)
        self.result_rows = -1
        self.result_size = -1
        self.secret_arn = secret_arn
        self.status = "STARTED"
        self.sub_statements: List[str] = []
        self.updated_at = now

    def __iter__(self) -> Iterator[Tuple[str, Any]]:
        yield "Id", self.id
        yield "ClusterIdentifier", self.cluster_identifier
        yield "CreatedAt", self.created_at
        yield "Database", self.database
        yield "DbUser", self.db_user
        yield "Duration", self.duration
        yield "HasResultSet", self.has_result_set
        yield "QueryParameters", self.query_parameters
        yield "QueryString", self.query_string
        yield "RedshiftPid", self.redshift_pid
        yield "RedshiftQueryId", self.redshift_query_id
        yield "ResultRows", self.result_rows
        yield "ResultSize", self.result_size
        yield "SecretArn", self.secret_arn
        yield "Status", self.status
        yield "SubStatements", self.sub_statements
        yield "UpdatedAt", self.updated_at


class StatementResult(Iterable[Tuple[str, Any]]):
    def __init__(
        self,
        column_metadata: List[Dict[str, Any]],
        records: List[List[Dict[str, Any]]],
        total_number_rows: int,
        next_token: Optional[str] = None,
    ):
        self.column_metadata = column_metadata
        self.records = records
        self.total_number_rows = total_number_rows
        self.next_token = next_token

    def __iter__(self) -> Iterator[Tuple[str, Any]]:
        yield "ColumnMetadata", self.column_metadata
        yield "Records", self.records
        yield "TotalNumberRows", self.total_number_rows
        yield "NextToken", self.next_token


class ColumnMetadata(Iterable[Tuple[str, Any]]):
    def __init__(
        self,
        column_default: Optional[str],
        is_case_sensitive: bool,
        is_signed: bool,
        name: str,
        nullable: int,
    ):
        self.column_default = column_default
        self.is_case_sensitive = is_case_sensitive
        self.is_signed = is_signed
        self.name = name
        self.nullable = nullable

    def __iter__(self) -> Iterator[Tuple[str, Any]]:
        yield "columnDefault", self.column_default
        yield "isCaseSensitive", self.is_case_sensitive
        yield "isSigned", self.is_signed
        yield "name", self.name
        yield "nullable", self.nullable


class Record(Iterable[Tuple[str, Any]]):
    def __init__(self, **kwargs: Any):
        self.kwargs = kwargs

    def __iter__(self) -> Iterator[Tuple[str, Any]]:
        if "long_value" in self.kwargs:
            yield "longValue", self.kwargs["long_value"]
        elif "string_value" in self.kwargs:
            yield "stringValue", self.kwargs["string_value"]


class RedshiftDataAPIServiceBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.statements: Dict[str, Statement] = {}

    @staticmethod
    def default_vpc_endpoint_service(
        service_region: str, zones: List[str]
    ) -> List[Dict[str, str]]:
        """Default VPC endpoint service."""
        return BaseBackend.default_vpc_endpoint_service_factory(
            service_region, zones, "redshift-data", policy_supported=False
        )

    def cancel_statement(self, statement_id: str) -> None:
        _validate_uuid(statement_id)

        try:
            # Statement exists
            statement = self.statements[statement_id]

            if statement.status != "STARTED":
                raise ValidationException(
                    "Could not cancel a query that is already in %s state with ID: %s"
                    % (statement.status, statement_id)
                )

            statement.status = "ABORTED"
            self.statements[statement_id] = statement
        except KeyError:
            # Statement does not exist.
            raise ResourceNotFoundException()

    def describe_statement(self, statement_id: str) -> Statement:
        _validate_uuid(statement_id)

        try:
            # Statement exists
            return self.statements[statement_id]
        except KeyError:
            # Statement does not exist.
            raise ResourceNotFoundException()

    def execute_statement(
        self,
        cluster_identifier: str,
        database: str,
        db_user: str,
        parameters: List[Dict[str, str]],
        secret_arn: str,
        sql: str,
    ) -> Statement:
        """
        Runs an SQL statement
        Validation of parameters is very limited because there is no redshift integration
        """
        statement = Statement(
            cluster_identifier=cluster_identifier,
            database=database,
            db_user=db_user,
            query_parameters=parameters,
            query_string=sql,
            secret_arn=secret_arn,
        )
        self.statements[statement.id] = statement
        return statement

    def get_statement_result(self, statement_id: str) -> StatementResult:
        """
        Return static statement result
        StatementResult is the result of the SQL query "sql" passed as parameter when calling "execute_statement"
        As such, it cannot be mocked
        """
        _validate_uuid(statement_id)

        if statement_id not in self.statements:
            raise ResourceNotFoundException()

        return StatementResult(
            [
                dict(ColumnMetadata(None, False, True, "Number", False)),
                dict(ColumnMetadata(None, True, False, "Street", False)),
                dict(ColumnMetadata(None, True, False, "City", False)),
            ],
            [
                [
                    dict(Record(long_value=10)),
                    dict(Record(string_value="Alpha st")),
                    dict(Record(string_value="Vancouver")),
                ],
                [
                    dict(Record(long_value=50)),
                    dict(Record(string_value="Beta st")),
                    dict(Record(string_value="Toronto")),
                ],
                [
                    dict(Record(long_value=100)),
                    dict(Record(string_value="Gamma av")),
                    dict(Record(string_value="Seattle")),
                ],
            ],
            3,
        )


def _validate_uuid(uuid: str) -> None:
    match = re.search(r"^[a-z0-9]{8}(-[a-z0-9]{4}){3}-[a-z0-9]{12}(:\d+)?$", uuid)
    if not match:
        raise ValidationException(
            "id must satisfy regex pattern: ^[a-z0-9]{8}(-[a-z0-9]{4}){3}-[a-z0-9]{12}(:\\d+)?$"
        )


# For unknown reasons I cannot use the service name "redshift-data" as I should
# It seems boto3 is unable to get the list of available regions for "redshift-data"
# See code here https://github.com/getmoto/moto/blob/master/moto/core/utils.py#L407
# sess.get_available_regions("redshift-data") returns an empty list
# Then I use the service redshift since they share the same regions
# See https://docs.aws.amazon.com/general/latest/gr/redshift-service.html
redshiftdata_backends = BackendDict(RedshiftDataAPIServiceBackend, "redshift")
