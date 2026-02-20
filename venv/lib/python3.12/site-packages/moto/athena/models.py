import time
from datetime import datetime
from typing import Any, Optional

from moto.athena.exceptions import InvalidArgumentException, QueryStillRunning
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.s3.models import s3_backends
from moto.s3.utils import bucket_and_name_from_url
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition


class TaggableResourceMixin:
    # This mixing was copied from Redshift when initially implementing
    # Athena. TBD if it's worth the overhead.

    def __init__(
        self,
        account_id: str,
        region_name: str,
        resource_name: str,
        tags: list[dict[str, str]],
    ):
        self.region = region_name
        self.resource_name = resource_name
        self.tags = tags or []
        self.arn = f"arn:{get_partition(region_name)}:athena:{region_name}:{account_id}:{resource_name}"

    def create_tags(self, tags: list[dict[str, str]]) -> list[dict[str, str]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys]
        self.tags.extend(tags)
        return self.tags

    def delete_tags(self, tag_keys: list[str]) -> list[dict[str, str]]:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]
        return self.tags


class WorkGroup(TaggableResourceMixin, BaseModel):
    resource_type = "workgroup"
    state = "ENABLED"

    def __init__(
        self,
        athena_backend: "AthenaBackend",
        name: str,
        configuration: dict[str, Any],
        description: str,
        tags: list[dict[str, str]],
    ):
        self.region_name = athena_backend.region_name
        super().__init__(
            athena_backend.account_id,
            self.region_name,
            f"workgroup/{name}",
            tags,
        )
        self.athena_backend = athena_backend
        self.name = name
        self.description = description
        self.configuration = configuration

        if "EnableMinimumEncryptionConfiguration" not in self.configuration:
            self.configuration["EnableMinimumEncryptionConfiguration"] = False
        if "EnforceWorkGroupConfiguration" not in self.configuration:
            self.configuration["EnforceWorkGroupConfiguration"] = True
        if "EngineVersion" not in self.configuration:
            self.configuration["EngineVersion"] = {
                "EffectiveEngineVersion": "Athena engine version 3",
                "SelectedEngineVersion": "AUTO",
            }
        if "PublishCloudWatchMetricsEnabled" not in self.configuration:
            self.configuration["PublishCloudWatchMetricsEnabled"] = False
        if "RequesterPaysEnabled" not in self.configuration:
            self.configuration["RequesterPaysEnabled"] = False


class DataCatalog(TaggableResourceMixin, BaseModel):
    def __init__(
        self,
        athena_backend: "AthenaBackend",
        name: str,
        catalog_type: str,
        description: str,
        parameters: str,
        tags: list[dict[str, str]],
    ):
        self.region_name = athena_backend.region_name
        super().__init__(
            athena_backend.account_id,
            self.region_name,
            f"datacatalog/{name}",
            tags,
        )
        self.athena_backend = athena_backend
        self.name = name
        self.type = catalog_type
        self.description = description
        self.parameters = parameters


class Execution(ManagedState):
    def __init__(
        self,
        query: str,
        context: str,
        config: dict[str, Any],
        workgroup: Optional[WorkGroup],
        execution_parameters: Optional[list[str]],
    ):
        ManagedState.__init__(
            self,
            model_name="athena::execution",
            transitions=[("QUEUED", "RUNNING"), ("RUNNING", "SUCCEEDED")],
        )
        self.id = str(mock_random.uuid4())
        self.query = query
        self.context = context
        self.config = config
        self.workgroup = workgroup
        self.execution_parameters = execution_parameters
        self.start_time = time.time()
        self.end_time = time.time()

        if self.config is not None and "OutputLocation" in self.config:
            if not self.config["OutputLocation"].endswith("/"):
                self.config["OutputLocation"] += "/"
            self.config["OutputLocation"] += f"{self.id}.csv"


class QueryResults(BaseModel):
    def __init__(self, rows: list[dict[str, Any]], column_info: list[dict[str, str]]):
        self.rows = rows
        self.column_info = column_info

    def to_dict(self) -> dict[str, Any]:
        return {
            "ResultSet": {
                "Rows": self.rows,
                "ResultSetMetadata": {"ColumnInfo": self.column_info},
            },
        }


class CapacityReservation(TaggableResourceMixin, BaseModel):
    def __init__(
        self,
        athena_backend: "AthenaBackend",
        name: str,
        target_dpus: int,
        tags: list[dict[str, str]],
    ):
        self.region_name = athena_backend.region_name
        super().__init__(
            athena_backend.account_id,
            self.region_name,
            f"capacity-reservation/{name}",
            tags,
        )
        self.athena_backend = athena_backend
        self.name = name
        self.target_dpus = target_dpus
        self.create_tags(tags)
        self.tags = tags


class NamedQuery(BaseModel):
    def __init__(
        self,
        name: str,
        description: str,
        database: str,
        query_string: str,
        workgroup: WorkGroup,
    ):
        self.id = str(mock_random.uuid4())
        self.name = name
        self.description = description
        self.database = database
        self.query_string = query_string
        self.workgroup = workgroup


class PreparedStatement(BaseModel):
    def __init__(
        self,
        statement_name: str,
        workgroup: WorkGroup,
        query_statement: str,
        description: str,
    ):
        self.statement_name = statement_name
        self.workgroup = workgroup
        self.query_statement = query_statement
        self.description = description
        self.last_modified_time = datetime.now()


class AthenaBackend(BaseBackend):
    PAGINATION_MODEL = {
        "list_named_queries": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 50,
            "unique_attribute": "id",
        }
    }

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.work_groups: dict[str, WorkGroup] = {}
        self.executions: dict[str, Execution] = {}
        self.named_queries: dict[str, NamedQuery] = {}
        self.capacity_reservations: dict[str, CapacityReservation] = {}
        self.data_catalogs: dict[str, DataCatalog] = {}
        self.query_results: dict[str, QueryResults] = {}
        self.query_results_queue: list[QueryResults] = []
        self.prepared_statements: dict[str, PreparedStatement] = {}
        self.tagger = TaggingService()

        # Initialise with the primary workgroup
        self.create_work_group(
            name="primary",
            description="",
            configuration={
                "ResultConfiguration": {},
                "EnforceWorkGroupConfiguration": False,
            },
            tags=[],
        )

    def create_work_group(
        self,
        name: str,
        configuration: dict[str, Any],
        description: str,
        tags: list[dict[str, str]],
    ) -> Optional[WorkGroup]:
        if name in self.work_groups:
            return None
        work_group = WorkGroup(self, name, configuration, description, tags)
        self.work_groups[name] = work_group
        self.tagger.tag_resource(work_group.arn, tags)
        return work_group

    def list_work_groups(self) -> list[dict[str, Any]]:
        return [
            {
                "Name": wg.name,
                "State": wg.state,
                "Description": wg.description,
                "CreationTime": time.time(),
            }
            for wg in self.work_groups.values()
        ]

    def get_work_group(self, name: str) -> Optional[dict[str, Any]]:
        if name not in self.work_groups:
            return None
        wg = self.work_groups[name]
        return {
            "Name": wg.name,
            "State": wg.state,
            "Configuration": wg.configuration,
            "Description": wg.description,
            "CreationTime": time.time(),
        }

    def delete_work_group(self, name: str) -> None:
        self.work_groups.pop(name, None)

    def start_query_execution(
        self,
        query: str,
        context: str,
        config: dict[str, Any],
        workgroup: str,
        execution_parameters: Optional[list[str]],
    ) -> str:
        execution = Execution(
            query=query,
            context=context,
            config=config,
            workgroup=self.work_groups.get(workgroup),
            execution_parameters=execution_parameters,
        )
        self.executions[execution.id] = execution

        self._store_predefined_query_results(execution.id)

        return execution.id

    def _store_predefined_query_results(self, exec_id: str) -> None:
        if exec_id not in self.query_results and self.query_results_queue:
            self.query_results[exec_id] = self.query_results_queue.pop(0)

            self._store_query_result_in_s3(exec_id)

    def get_query_execution(self, exec_id: str) -> Execution:
        execution = self.executions[exec_id]
        execution.advance()
        return execution

    def list_query_executions(self, workgroup: Optional[str]) -> dict[str, Execution]:
        # Note: We do not advance the execution status here, only in `get_query_execution`
        # This method simply returns the QueryExecutionIds to the user
        # They will always have to call `get_query_execution` to get the status
        if workgroup is not None:
            return {
                exec_id: execution
                for exec_id, execution in self.executions.items()
                if execution.workgroup and execution.workgroup.name == workgroup
            }
        return self.executions

    def get_query_results(self, exec_id: str) -> QueryResults:
        """
        Queries are not executed by Moto, so this call will always return 0 rows by default.

        You can use a dedicated API to override this, by configuring a queue of expected results.

        A request to `get_query_results` will take the first result from that queue, and assign it to the provided QueryExecutionId. Subsequent requests using the same QueryExecutionId will return the same result. Other requests using a different QueryExecutionId will take the next result from the queue, or return an empty result if the queue is empty.

        Configuring this queue by making an HTTP request to `/moto-api/static/athena/query-results`. An example invocation looks like this:

        .. sourcecode:: python

            expected_results = {
                "account_id": "123456789012",  # This is the default - can be omitted
                "region": "us-east-1",  # This is the default - can be omitted
                "results": [
                    {
                        "rows": [{"Data": [{"VarCharValue": "1"}]}],
                        "column_info": [{
                            "CatalogName": "string",
                            "SchemaName": "string",
                            "TableName": "string",
                            "Name": "string",
                            "Label": "string",
                            "Type": "string",
                            "Precision": 123,
                            "Scale": 123,
                            "Nullable": "NOT_NULL",
                            "CaseSensitive": True,
                        }],
                    },
                    # other results as required
                ],
            }
            resp = requests.post(
                "http://motoapi.amazonaws.com/moto-api/static/athena/query-results",
                json=expected_results,
            )
            assert resp.status_code == 201

            client = boto3.client("athena", region_name="us-east-1")
            details = client.get_query_execution(QueryExecutionId="any_id")["QueryExecution"]

        .. note:: The exact QueryExecutionId is not relevant here, but will likely be whatever value is returned by start_query_execution

        Query results will also be stored in the S3 output location (in CSV format).

        """
        if (exctn := self.executions.get(exec_id)) and exctn.status != "SUCCEEDED":
            raise QueryStillRunning(current_status=exctn.status)

        self._store_predefined_query_results(exec_id)

        results = (
            self.query_results[exec_id]
            if exec_id in self.query_results
            else QueryResults(rows=[], column_info=[])
        )
        return results

    def _store_query_result_in_s3(self, exec_id: str) -> None:
        try:
            output_location = self.executions[exec_id].config["OutputLocation"]
            bucket, key = bucket_and_name_from_url(output_location)

            query_result = ""
            for row in self.query_results[exec_id].rows:
                query_result += ",".join(
                    [
                        f'"{r["VarCharValue"]}"' if "VarCharValue" in r else ""
                        for r in row["Data"]
                    ]
                )
                query_result += "\n"

            s3_backends[self.account_id][self.partition].put_object(
                bucket_name=bucket,  # type: ignore
                key_name=key,  # type: ignore
                value=query_result.encode("utf-8"),
            )
        except:  # noqa
            # Execution may not exist
            # OutputLocation may not exist
            pass

    def stop_query_execution(self, exec_id: str) -> None:
        execution = self.executions[exec_id]
        execution.status = "CANCELLED"

    def create_capacity_reservation(
        self,
        name: str,
        target_dpus: int,
        tags: list[dict[str, str]],
    ) -> None:
        cr = CapacityReservation(self, name, target_dpus, tags)
        self.capacity_reservations[cr.name] = cr
        self.tagger.tag_resource(cr.arn, tags)
        return None

    def get_capacity_reservation(self, name: str) -> Optional[CapacityReservation]:
        return self.capacity_reservations.get(name)

    def list_capacity_reservations(self) -> list[dict[str, Any]]:
        return [
            {"Name": cr.name, "TargetDpus": cr.target_dpus, "CreationTime": time.time()}
            for cr in self.capacity_reservations.values()
        ]

    def update_capacity_reservation(self, name: str, target_dpus: int) -> None:
        if name not in self.capacity_reservations:
            raise InvalidArgumentException("Capacity Reservation does not exist")

        self.capacity_reservations[name].target_dpus = target_dpus

    def create_named_query(
        self,
        name: str,
        description: str,
        database: str,
        query_string: str,
        workgroup: str,
    ) -> str:
        nq = NamedQuery(
            name=name,
            description=description,
            database=database,
            query_string=query_string,
            workgroup=self.work_groups[workgroup],
        )
        self.named_queries[nq.id] = nq
        return nq.id

    def get_named_query(self, query_id: str) -> Optional[NamedQuery]:
        return self.named_queries[query_id] if query_id in self.named_queries else None

    def list_data_catalogs(self) -> list[dict[str, str]]:
        return [
            {"CatalogName": dc.name, "Type": dc.type}
            for dc in self.data_catalogs.values()
        ]

    def get_data_catalog(self, name: str) -> Optional[dict[str, str]]:
        if name not in self.data_catalogs:
            return None
        dc = self.data_catalogs[name]
        return {
            "Name": dc.name,
            "Description": dc.description,
            "Type": dc.type,
            "Parameters": dc.parameters,
        }

    def create_data_catalog(
        self,
        name: str,
        catalog_type: str,
        description: str,
        parameters: str,
        tags: list[dict[str, str]],
    ) -> Optional[DataCatalog]:
        if name in self.data_catalogs:
            return None
        data_catalog = DataCatalog(
            self, name, catalog_type, description, parameters, tags
        )
        self.data_catalogs[name] = data_catalog
        self.tagger.tag_resource(data_catalog.arn, tags)
        return data_catalog

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_named_queries(self, work_group: str) -> list[str]:
        named_query_ids = [
            q.id for q in self.named_queries.values() if q.workgroup.name == work_group
        ]
        return named_query_ids

    def create_prepared_statement(
        self,
        statement_name: str,
        workgroup: WorkGroup,
        query_statement: str,
        description: str,
    ) -> None:
        ps = PreparedStatement(
            statement_name=statement_name,
            workgroup=workgroup,
            query_statement=query_statement,
            description=description,
        )
        self.prepared_statements[ps.statement_name] = ps
        return None

    def get_prepared_statement(
        self, statement_name: str, work_group: WorkGroup
    ) -> Optional[PreparedStatement]:
        if statement_name in self.prepared_statements:
            ps = self.prepared_statements[statement_name]
            if ps.workgroup == work_group:
                return ps
        return None

    def get_query_runtime_statistics(
        self, query_execution_id: str
    ) -> Optional[Execution]:
        if query_execution_id in self.executions:
            return self.executions[query_execution_id]
        return None

    def list_tags_for_resource(self, resource_arn: str) -> Optional[dict[str, Any]]:
        if self.tagger.has_tags(resource_arn):
            return self.tagger.list_tags_for_resource(resource_arn)
        return None


athena_backends = BackendDict(AthenaBackend, "athena")
