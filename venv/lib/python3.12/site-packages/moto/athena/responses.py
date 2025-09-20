import json
from typing import Dict, Tuple, Union

from moto.core.responses import BaseResponse

from .models import AthenaBackend, athena_backends


class AthenaResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="athena")

    @property
    def athena_backend(self) -> AthenaBackend:
        return athena_backends[self.current_account][self.region]

    def create_work_group(self) -> Union[Tuple[str, Dict[str, int]], str]:
        name = self._get_param("Name")
        description = self._get_param("Description")
        configuration = self._get_param("Configuration")
        tags = self._get_param("Tags")
        work_group = self.athena_backend.create_work_group(
            name, configuration, description, tags
        )
        if not work_group:
            return self.error("WorkGroup already exists", 400)
        return json.dumps(
            {
                "CreateWorkGroupResponse": {
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    }
                }
            }
        )

    def list_work_groups(self) -> str:
        return json.dumps({"WorkGroups": self.athena_backend.list_work_groups()})

    def get_work_group(self) -> str:
        name = self._get_param("WorkGroup")
        return json.dumps({"WorkGroup": self.athena_backend.get_work_group(name)})

    def delete_work_group(self) -> str:
        name = self._get_param("WorkGroup")
        self.athena_backend.delete_work_group(name)
        return "{}"

    def start_query_execution(self) -> Union[Tuple[str, Dict[str, int]], str]:
        query = self._get_param("QueryString")
        context = self._get_param("QueryExecutionContext")
        config = self._get_param("ResultConfiguration")
        workgroup = self._get_param("WorkGroup")
        execution_parameters = self._get_param("ExecutionParameters")
        if workgroup and not self.athena_backend.get_work_group(workgroup):
            return self.error("WorkGroup does not exist", 400)
        q_exec_id = self.athena_backend.start_query_execution(
            query=query,
            context=context,
            config=config,
            workgroup=workgroup,
            execution_parameters=execution_parameters,
        )
        return json.dumps({"QueryExecutionId": q_exec_id})

    def get_query_execution(self) -> str:
        exec_id = self._get_param("QueryExecutionId")
        execution = self.athena_backend.get_query_execution(exec_id)
        ddl_commands = ("ALTER", "CREATE", "DESCRIBE", "DROP", "MSCK", "SHOW")
        statement_type = "DML"
        if execution.query.upper().startswith(ddl_commands):
            statement_type = "DDL"
        result = {
            "QueryExecution": {
                "QueryExecutionId": exec_id,
                "Query": execution.query,
                "StatementType": statement_type,
                "ResultConfiguration": execution.config,
                "ResultReuseConfiguration": {
                    "ResultReuseByAgeConfiguration": {"Enabled": False}
                },
                "QueryExecutionContext": execution.context,
                "Status": {
                    "State": execution.status,
                    "SubmissionDateTime": execution.start_time,
                    "CompletionDateTime": execution.end_time,
                },
                "Statistics": {
                    "EngineExecutionTimeInMillis": 0,
                    "DataScannedInBytes": 0,
                    "TotalExecutionTimeInMillis": 0,
                    "QueryQueueTimeInMillis": 0,
                    "ServicePreProcessingTimeInMillis": 0,
                    "QueryPlanningTimeInMillis": 0,
                    "ServiceProcessingTimeInMillis": 0,
                    "ResultReuseInformation": {"ReusedPreviousResult": False},
                },
                "WorkGroup": execution.workgroup.name if execution.workgroup else None,
            }
        }
        if execution.execution_parameters is not None:
            result["QueryExecution"]["ExecutionParameters"] = (
                execution.execution_parameters
            )
        return json.dumps(result)

    def create_capacity_reservation(self) -> Union[Tuple[str, Dict[str, int]], str]:
        name = self._get_param("Name")
        target_dpus = self._get_param("TargetDpus")
        tags = self._get_param("Tags")
        self.athena_backend.create_capacity_reservation(name, target_dpus, tags)
        return json.dumps(dict())

    def get_capacity_reservation(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        name = self._get_param("Name")
        capacity_reservation = self.athena_backend.get_capacity_reservation(name)
        if not capacity_reservation:
            return self.error("Capacity reservation does not exist", 400)
        return json.dumps(
            {
                "CapacityReservation": {
                    "Name": capacity_reservation.name,
                    "TargetDpus": capacity_reservation.target_dpus,
                    "Tags": capacity_reservation.tags,
                }
            }
        )

    def get_query_results(self) -> str:
        exec_id = self._get_param("QueryExecutionId")
        result = self.athena_backend.get_query_results(exec_id)
        return json.dumps(result.to_dict())

    def list_query_executions(self) -> str:
        workgroup = self._get_param("WorkGroup")
        executions = self.athena_backend.list_query_executions(workgroup)
        return json.dumps({"QueryExecutionIds": [i for i in executions.keys()]})

    def stop_query_execution(self) -> str:
        exec_id = self._get_param("QueryExecutionId")
        self.athena_backend.stop_query_execution(exec_id)
        return json.dumps({})

    def error(self, msg: str, status: int) -> Tuple[str, Dict[str, int]]:
        return (
            json.dumps({"__type": "InvalidRequestException", "Message": msg}),
            dict(status=status),
        )

    def create_named_query(self) -> Union[Tuple[str, Dict[str, int]], str]:
        name = self._get_param("Name")
        description = self._get_param("Description")
        database = self._get_param("Database")
        query_string = self._get_param("QueryString")
        workgroup = self._get_param("WorkGroup") or "primary"
        if not self.athena_backend.get_work_group(workgroup):
            return self.error("WorkGroup does not exist", 400)
        query_id = self.athena_backend.create_named_query(
            name, description, database, query_string, workgroup
        )
        return json.dumps({"NamedQueryId": query_id})

    def get_named_query(self) -> str:
        query_id = self._get_param("NamedQueryId")
        nq = self.athena_backend.get_named_query(query_id)
        return json.dumps(
            {
                "NamedQuery": {
                    "Name": nq.name,  # type: ignore[union-attr]
                    "Description": nq.description,  # type: ignore[union-attr]
                    "Database": nq.database,  # type: ignore[union-attr]
                    "QueryString": nq.query_string,  # type: ignore[union-attr]
                    "NamedQueryId": nq.id,  # type: ignore[union-attr]
                    "WorkGroup": nq.workgroup.name,  # type: ignore[union-attr]
                }
            }
        )

    def list_data_catalogs(self) -> str:
        return json.dumps(
            {"DataCatalogsSummary": self.athena_backend.list_data_catalogs()}
        )

    def list_tags_for_resource(self) -> Union[Tuple[str, Dict[str, int]], str]:
        resource_arn = self._get_param("ResourceARN")
        tags = self.athena_backend.list_tags_for_resource(resource_arn)
        if not tags:
            return self.error(f"Athena Resource, {resource_arn} Does Not Exist", 400)
        return json.dumps(tags)

    def get_data_catalog(self) -> str:
        name = self._get_param("Name")
        return json.dumps({"DataCatalog": self.athena_backend.get_data_catalog(name)})

    def create_data_catalog(self) -> Union[Tuple[str, Dict[str, int]], str]:
        name = self._get_param("Name")
        catalog_type = self._get_param("Type")
        description = self._get_param("Description")
        parameters = self._get_param("Parameters")
        tags = self._get_param("Tags")
        data_catalog = self.athena_backend.create_data_catalog(
            name, catalog_type, description, parameters, tags
        )
        if not data_catalog:
            return self.error("DataCatalog already exists", 400)
        return json.dumps(
            {
                "CreateDataCatalogResponse": {
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    }
                }
            }
        )

    def list_named_queries(self) -> str:
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")
        work_group = self._get_param("WorkGroup") or "primary"
        named_query_ids, next_token = self.athena_backend.list_named_queries(
            next_token=next_token, max_results=max_results, work_group=work_group
        )
        return json.dumps({"NamedQueryIds": named_query_ids, "NextToken": next_token})

    def create_prepared_statement(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        statement_name = self._get_param("StatementName")
        work_group = self._get_param("WorkGroup")
        query_statement = self._get_param("QueryStatement")
        description = self._get_param("Description")
        if not self.athena_backend.get_work_group(work_group):
            return self.error("WorkGroup does not exist", 400)
        self.athena_backend.create_prepared_statement(
            statement_name=statement_name,
            workgroup=work_group,
            query_statement=query_statement,
            description=description,
        )
        return json.dumps(dict())

    def get_prepared_statement(self) -> str:
        statement_name = self._get_param("StatementName")
        work_group = self._get_param("WorkGroup")
        ps = self.athena_backend.get_prepared_statement(
            statement_name=statement_name,
            work_group=work_group,
        )
        return json.dumps(
            {
                "PreparedStatement": {
                    "StatementName": ps.statement_name,  # type: ignore[union-attr]
                    "QueryStatement": ps.query_statement,  # type: ignore[union-attr]
                    "WorkGroupName": ps.workgroup,  # type: ignore[union-attr]
                    "Description": ps.description,  # type: ignore[union-attr]
                    # "LastModifiedTime": ps.last_modified_time,  # type: ignore[union-attr]
                }
            }
        )

    def get_query_runtime_statistics(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        query_execution_id = self._get_param("QueryExecutionId")

        ps = self.athena_backend.get_query_runtime_statistics(
            query_execution_id=query_execution_id
        )

        if ps is None:
            return self.error(f"QueryExecution {query_execution_id} was not found", 400)

        return json.dumps(
            {
                "QueryRuntimeStatistics": {
                    "OutputStage": {
                        "ExecutionTime": 100,
                        "InputBytes": 0,
                        "InputRows": 0,
                        "OutputBytes": 1,
                        "OutputRows": 1,
                        "StageId": 1,
                        "State": ps.status,
                    },
                    "Rows": {
                        "InputBytes": 0,
                        "InputRows": 0,
                        "OutputBytes": 2,
                        "OutputRows": 2,
                    },
                    "Timeline": {
                        "EngineExecutionTimeInMillis": 0,
                        "QueryPlanningTimeInMillis": 0,
                        "QueryQueueTimeInMillis": 0,
                        "ServicePreProcessingTimeInMillis": 0,
                        "ServiceProcessingTimeInMillis": 0,
                        "TotalExecutionTimeInMillis": 0,
                    },
                }
            }
        )
