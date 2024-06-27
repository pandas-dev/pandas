import time
from datetime import datetime
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.utilities.paginator import paginate
from moto.utilities.utils import get_partition


class TaggableResourceMixin:
    # This mixing was copied from Redshift when initially implementing
    # Athena. TBD if it's worth the overhead.

    def __init__(
        self,
        account_id: str,
        region_name: str,
        resource_name: str,
        tags: List[Dict[str, str]],
    ):
        self.region = region_name
        self.resource_name = resource_name
        self.tags = tags or []
        self.arn = f"arn:{get_partition(region_name)}:athena:{region_name}:{account_id}:{resource_name}"

    def create_tags(self, tags: List[Dict[str, str]]) -> List[Dict[str, str]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys]
        self.tags.extend(tags)
        return self.tags

    def delete_tags(self, tag_keys: List[str]) -> List[Dict[str, str]]:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]
        return self.tags


class WorkGroup(TaggableResourceMixin, BaseModel):
    resource_type = "workgroup"
    state = "ENABLED"

    def __init__(
        self,
        athena_backend: "AthenaBackend",
        name: str,
        configuration: Dict[str, Any],
        description: str,
        tags: List[Dict[str, str]],
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


class DataCatalog(TaggableResourceMixin, BaseModel):
    def __init__(
        self,
        athena_backend: "AthenaBackend",
        name: str,
        catalog_type: str,
        description: str,
        parameters: str,
        tags: List[Dict[str, str]],
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


class Execution(BaseModel):
    def __init__(
        self,
        query: str,
        context: str,
        config: Dict[str, Any],
        workgroup: WorkGroup,
        execution_parameters: Optional[List[str]],
    ):
        self.id = str(mock_random.uuid4())
        self.query = query
        self.context = context
        self.config = config
        self.workgroup = workgroup
        self.execution_parameters = execution_parameters
        self.start_time = time.time()
        self.status = "SUCCEEDED"

        if self.config is not None and "OutputLocation" in self.config:
            if not self.config["OutputLocation"].endswith("/"):
                self.config["OutputLocation"] += "/"
            self.config["OutputLocation"] += f"{self.id}.csv"


class QueryResults(BaseModel):
    def __init__(self, rows: List[Dict[str, Any]], column_info: List[Dict[str, str]]):
        self.rows = rows
        self.column_info = column_info

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ResultSet": {
                "Rows": self.rows,
                "ResultSetMetadata": {"ColumnInfo": self.column_info},
            },
        }


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
        self.work_groups: Dict[str, WorkGroup] = {}
        self.executions: Dict[str, Execution] = {}
        self.named_queries: Dict[str, NamedQuery] = {}
        self.data_catalogs: Dict[str, DataCatalog] = {}
        self.query_results: Dict[str, QueryResults] = {}
        self.query_results_queue: List[QueryResults] = []
        self.prepared_statements: Dict[str, PreparedStatement] = {}

        # Initialise with the primary workgroup
        self.create_work_group(
            name="primary", description="", configuration=dict(), tags=[]
        )

    def create_work_group(
        self,
        name: str,
        configuration: Dict[str, Any],
        description: str,
        tags: List[Dict[str, str]],
    ) -> Optional[WorkGroup]:
        if name in self.work_groups:
            return None
        work_group = WorkGroup(self, name, configuration, description, tags)
        self.work_groups[name] = work_group
        return work_group

    def list_work_groups(self) -> List[Dict[str, Any]]:
        return [
            {
                "Name": wg.name,
                "State": wg.state,
                "Description": wg.description,
                "CreationTime": time.time(),
            }
            for wg in self.work_groups.values()
        ]

    def get_work_group(self, name: str) -> Optional[Dict[str, Any]]:
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

    def start_query_execution(
        self,
        query: str,
        context: str,
        config: Dict[str, Any],
        workgroup: WorkGroup,
        execution_parameters: Optional[List[str]],
    ) -> str:
        execution = Execution(
            query=query,
            context=context,
            config=config,
            workgroup=workgroup,
            execution_parameters=execution_parameters,
        )
        self.executions[execution.id] = execution
        return execution.id

    def get_query_execution(self, exec_id: str) -> Execution:
        return self.executions[exec_id]

    def list_query_executions(self) -> Dict[str, Execution]:
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

        """
        if exec_id not in self.query_results and self.query_results_queue:
            self.query_results[exec_id] = self.query_results_queue.pop(0)
        results = (
            self.query_results[exec_id]
            if exec_id in self.query_results
            else QueryResults(rows=[], column_info=[])
        )
        return results

    def stop_query_execution(self, exec_id: str) -> None:
        execution = self.executions[exec_id]
        execution.status = "CANCELLED"

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

    def list_data_catalogs(self) -> List[Dict[str, str]]:
        return [
            {"CatalogName": dc.name, "Type": dc.type}
            for dc in self.data_catalogs.values()
        ]

    def get_data_catalog(self, name: str) -> Optional[Dict[str, str]]:
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
        tags: List[Dict[str, str]],
    ) -> Optional[DataCatalog]:
        if name in self.data_catalogs:
            return None
        data_catalog = DataCatalog(
            self, name, catalog_type, description, parameters, tags
        )
        self.data_catalogs[name] = data_catalog
        return data_catalog

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_named_queries(self, work_group: str) -> List[str]:
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


athena_backends = BackendDict(AthenaBackend, "athena")
