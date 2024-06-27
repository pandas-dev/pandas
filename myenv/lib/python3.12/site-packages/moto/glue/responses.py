import json
from typing import Any, Dict, List

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import (
    FakeCrawler,
    FakeJob,
    FakeSession,
    FakeTrigger,
    GlueBackend,
    glue_backends,
)


class GlueResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="glue")

    @property
    def glue_backend(self) -> GlueBackend:
        return glue_backends[self.current_account][self.region]

    @property
    def parameters(self) -> Dict[str, Any]:  # type: ignore[misc]
        return json.loads(self.body)

    def create_database(self) -> str:
        database_input = self.parameters.get("DatabaseInput")
        database_name = database_input.get("Name")  # type: ignore
        if "CatalogId" in self.parameters:
            database_input["CatalogId"] = self.parameters.get("CatalogId")  # type: ignore
        self.glue_backend.create_database(
            database_name,
            database_input,  # type: ignore[arg-type]
            self.parameters.get("Tags"),  # type: ignore[arg-type]
        )
        return ""

    def get_database(self) -> str:
        database_name = self.parameters.get("Name")
        database = self.glue_backend.get_database(database_name)  # type: ignore[arg-type]
        return json.dumps({"Database": database.as_dict()})

    def get_databases(self) -> str:
        database_list = self.glue_backend.get_databases()
        return json.dumps(
            {"DatabaseList": [database.as_dict() for database in database_list]}
        )

    def update_database(self) -> str:
        database_input = self.parameters.get("DatabaseInput")
        database_name = self.parameters.get("Name")
        if "CatalogId" in self.parameters:
            database_input["CatalogId"] = self.parameters.get("CatalogId")  # type: ignore
        self.glue_backend.update_database(database_name, database_input)  # type: ignore[arg-type]
        return ""

    def delete_database(self) -> str:
        name = self.parameters.get("Name")
        self.glue_backend.delete_database(name)  # type: ignore[arg-type]
        return json.dumps({})

    def create_table(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_input = self.parameters.get("TableInput")
        table_name = table_input.get("Name")  # type: ignore
        self.glue_backend.create_table(database_name, table_name, table_input)  # type: ignore[arg-type]
        return ""

    def get_table(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("Name")
        table = self.glue_backend.get_table(database_name, table_name)  # type: ignore[arg-type]

        return json.dumps({"Table": table.as_dict()})

    def update_table(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_input = self.parameters.get("TableInput")
        table_name = table_input.get("Name")  # type: ignore
        self.glue_backend.update_table(database_name, table_name, table_input)  # type: ignore[arg-type]
        return ""

    def get_table_versions(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("TableName")
        versions = self.glue_backend.get_table_versions(database_name, table_name)  # type: ignore[arg-type]
        return json.dumps(
            {
                "TableVersions": [
                    {"Table": data, "VersionId": version}
                    for version, data in versions.items()
                ]
            }
        )

    def get_table_version(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("TableName")
        ver_id = self.parameters.get("VersionId")
        return self.glue_backend.get_table_version(database_name, table_name, ver_id)  # type: ignore[arg-type]

    def delete_table_version(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("TableName")
        version_id = self.parameters.get("VersionId")
        self.glue_backend.delete_table_version(database_name, table_name, version_id)  # type: ignore[arg-type]
        return "{}"

    def get_tables(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        expression = self.parameters.get("Expression")
        tables = self.glue_backend.get_tables(database_name, expression)  # type: ignore[arg-type]
        return json.dumps({"TableList": [table.as_dict() for table in tables]})

    def delete_table(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("Name")
        self.glue_backend.delete_table(database_name, table_name)  # type: ignore[arg-type]
        return "{}"

    def batch_delete_table(self) -> str:
        database_name = self.parameters.get("DatabaseName")

        tables = self.parameters.get("TablesToDelete")
        errors = self.glue_backend.batch_delete_table(database_name, tables)  # type: ignore[arg-type]

        out = {}
        if errors:
            out["Errors"] = errors

        return json.dumps(out)

    def get_partitions(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("TableName")
        expression = self.parameters.get("Expression")
        partitions = self.glue_backend.get_partitions(
            database_name,  # type: ignore[arg-type]
            table_name,  # type: ignore[arg-type]
            expression,  # type: ignore[arg-type]
        )

        return json.dumps({"Partitions": [p.as_dict() for p in partitions]})

    def get_partition(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("TableName")
        values = self.parameters.get("PartitionValues")

        p = self.glue_backend.get_partition(database_name, table_name, values)  # type: ignore[arg-type]

        return json.dumps({"Partition": p.as_dict()})

    def batch_get_partition(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("TableName")
        partitions_to_get = self.parameters.get("PartitionsToGet")

        partitions = self.glue_backend.batch_get_partition(
            database_name,  # type: ignore[arg-type]
            table_name,  # type: ignore[arg-type]
            partitions_to_get,  # type: ignore[arg-type]
        )

        return json.dumps({"Partitions": partitions})

    def create_partition(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("TableName")
        part_input = self.parameters.get("PartitionInput")

        self.glue_backend.create_partition(database_name, table_name, part_input)  # type: ignore[arg-type]
        return ""

    def batch_create_partition(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("TableName")
        partition_input = self.parameters.get("PartitionInputList")
        errors_output = self.glue_backend.batch_create_partition(
            database_name,  # type: ignore[arg-type]
            table_name,  # type: ignore[arg-type]
            partition_input,  # type: ignore[arg-type]
        )

        out = {}
        if errors_output:
            out["Errors"] = errors_output

        return json.dumps(out)

    def update_partition(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("TableName")
        part_input = self.parameters.get("PartitionInput")
        part_to_update = self.parameters.get("PartitionValueList")

        self.glue_backend.update_partition(
            database_name,  # type: ignore[arg-type]
            table_name,  # type: ignore[arg-type]
            part_input,  # type: ignore[arg-type]
            part_to_update,  # type: ignore[arg-type]
        )
        return ""

    def batch_update_partition(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("TableName")
        entries = self.parameters.get("Entries")

        errors_output = self.glue_backend.batch_update_partition(
            database_name,  # type: ignore[arg-type]
            table_name,  # type: ignore[arg-type]
            entries,  # type: ignore[arg-type]
        )

        out = {}
        if errors_output:
            out["Errors"] = errors_output

        return json.dumps(out)

    def delete_partition(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("TableName")
        part_to_delete = self.parameters.get("PartitionValues")

        self.glue_backend.delete_partition(database_name, table_name, part_to_delete)  # type: ignore[arg-type]
        return ""

    def batch_delete_partition(self) -> str:
        database_name = self.parameters.get("DatabaseName")
        table_name = self.parameters.get("TableName")
        parts = self.parameters.get("PartitionsToDelete")

        errors_output = self.glue_backend.batch_delete_partition(
            database_name,  # type: ignore[arg-type]
            table_name,  # type: ignore[arg-type]
            parts,  # type: ignore[arg-type]
        )

        out = {}
        if errors_output:
            out["Errors"] = errors_output

        return json.dumps(out)

    def create_crawler(self) -> str:
        self.glue_backend.create_crawler(
            name=self.parameters.get("Name"),  # type: ignore[arg-type]
            role=self.parameters.get("Role"),  # type: ignore[arg-type]
            database_name=self.parameters.get("DatabaseName"),  # type: ignore[arg-type]
            description=self.parameters.get("Description"),  # type: ignore[arg-type]
            targets=self.parameters.get("Targets"),  # type: ignore[arg-type]
            schedule=self.parameters.get("Schedule"),  # type: ignore[arg-type]
            classifiers=self.parameters.get("Classifiers"),  # type: ignore[arg-type]
            table_prefix=self.parameters.get("TablePrefix"),  # type: ignore[arg-type]
            schema_change_policy=self.parameters.get("SchemaChangePolicy"),  # type: ignore[arg-type]
            recrawl_policy=self.parameters.get("RecrawlPolicy"),  # type: ignore[arg-type]
            lineage_configuration=self.parameters.get("LineageConfiguration"),  # type: ignore[arg-type]
            configuration=self.parameters.get("Configuration"),  # type: ignore[arg-type]
            crawler_security_configuration=self.parameters.get(  # type: ignore[arg-type]
                "CrawlerSecurityConfiguration"
            ),
            tags=self.parameters.get("Tags"),  # type: ignore[arg-type]
        )
        return ""

    def get_crawler(self) -> str:
        name = self.parameters.get("Name")
        crawler = self.glue_backend.get_crawler(name)  # type: ignore[arg-type]
        return json.dumps({"Crawler": crawler.as_dict()})

    def get_crawlers(self) -> str:
        crawlers = self.glue_backend.get_crawlers()
        return json.dumps({"Crawlers": [crawler.as_dict() for crawler in crawlers]})

    def list_crawlers(self) -> str:
        next_token = self._get_param("NextToken")
        max_results = self._get_int_param("MaxResults")
        tags = self._get_param("Tags")
        crawlers, next_token = self.glue_backend.list_crawlers(
            next_token=next_token, max_results=max_results
        )
        filtered_crawler_names = self.filter_crawlers_by_tags(crawlers, tags)
        return json.dumps(
            dict(
                CrawlerNames=[crawler_name for crawler_name in filtered_crawler_names],
                NextToken=next_token,
            )
        )

    def filter_crawlers_by_tags(
        self, crawlers: List[FakeCrawler], tags: Dict[str, str]
    ) -> List[str]:
        if not tags:
            return [crawler.get_name() for crawler in crawlers]
        return [
            crawler.get_name()
            for crawler in crawlers
            if self.is_tags_match(crawler.arn, tags)
        ]

    def start_crawler(self) -> str:
        name = self.parameters.get("Name")
        self.glue_backend.start_crawler(name)  # type: ignore[arg-type]
        return ""

    def stop_crawler(self) -> str:
        name = self.parameters.get("Name")
        self.glue_backend.stop_crawler(name)  # type: ignore[arg-type]
        return ""

    def delete_crawler(self) -> str:
        name = self.parameters.get("Name")
        self.glue_backend.delete_crawler(name)  # type: ignore[arg-type]
        return ""

    def create_job(self) -> str:
        name = self._get_param("Name")
        description = self._get_param("Description")
        log_uri = self._get_param("LogUri")
        role = self._get_param("Role")
        execution_property = self._get_param("ExecutionProperty")
        command = self._get_param("Command")
        default_arguments = self._get_param("DefaultArguments")
        non_overridable_arguments = self._get_param("NonOverridableArguments")
        connections = self._get_param("Connections")
        max_retries = self._get_int_param("MaxRetries")
        allocated_capacity = self._get_int_param("AllocatedCapacity")
        timeout = self._get_int_param("Timeout")
        max_capacity = self._get_param("MaxCapacity")
        security_configuration = self._get_param("SecurityConfiguration")
        tags = self._get_param("Tags")
        notification_property = self._get_param("NotificationProperty")
        glue_version = self._get_param("GlueVersion")
        number_of_workers = self._get_int_param("NumberOfWorkers")
        worker_type = self._get_param("WorkerType")
        code_gen_configuration_nodes = self._get_param("CodeGenConfigurationNodes")
        execution_class = self._get_param("ExecutionClass")
        source_control_details = self._get_param("SourceControlDetails")
        self.glue_backend.create_job(
            name=name,
            description=description,
            log_uri=log_uri,
            role=role,
            execution_property=execution_property,
            command=command,
            default_arguments=default_arguments,
            non_overridable_arguments=non_overridable_arguments,
            connections=connections,
            max_retries=max_retries,
            allocated_capacity=allocated_capacity,
            timeout=timeout,
            max_capacity=max_capacity,
            security_configuration=security_configuration,
            tags=tags,
            notification_property=notification_property,
            glue_version=glue_version,
            number_of_workers=number_of_workers,
            worker_type=worker_type,
            code_gen_configuration_nodes=code_gen_configuration_nodes,
            execution_class=execution_class,
            source_control_details=source_control_details,
        )
        return json.dumps(dict(Name=name))

    def get_job(self) -> str:
        name = self.parameters.get("JobName")
        job = self.glue_backend.get_job(name)  # type: ignore[arg-type]
        return json.dumps({"Job": job.as_dict()})

    def get_jobs(self) -> str:
        next_token = self._get_param("NextToken")
        max_results = self._get_int_param("MaxResults")
        jobs, next_token = self.glue_backend.get_jobs(
            next_token=next_token, max_results=max_results
        )
        return json.dumps(
            dict(
                Jobs=[job.as_dict() for job in jobs],
                NextToken=next_token,
            )
        )

    def start_job_run(self) -> str:
        name = self.parameters.get("JobName")
        job_run_id = self.glue_backend.start_job_run(name)  # type: ignore[arg-type]
        return json.dumps(dict(JobRunId=job_run_id))

    def get_job_run(self) -> str:
        name = self.parameters.get("JobName")
        run_id = self.parameters.get("RunId")
        job_run = self.glue_backend.get_job_run(name, run_id)  # type: ignore[arg-type]
        return json.dumps({"JobRun": job_run.as_dict()})

    def list_jobs(self) -> str:
        next_token = self._get_param("NextToken")
        max_results = self._get_int_param("MaxResults")
        tags = self._get_param("Tags")
        jobs, next_token = self.glue_backend.list_jobs(
            next_token=next_token, max_results=max_results
        )
        filtered_job_names = self.filter_jobs_by_tags(jobs, tags)
        return json.dumps(
            dict(
                JobNames=[job_name for job_name in filtered_job_names],
                NextToken=next_token,
            )
        )

    def delete_job(self) -> str:
        name = self.parameters.get("JobName")
        self.glue_backend.delete_job(name)  # type: ignore[arg-type]
        return json.dumps({"JobName": name})

    def get_tags(self) -> TYPE_RESPONSE:
        resource_arn = self.parameters.get("ResourceArn")
        tags = self.glue_backend.get_tags(resource_arn)  # type: ignore[arg-type]
        return 200, {}, json.dumps({"Tags": tags})

    def tag_resource(self) -> TYPE_RESPONSE:
        resource_arn = self.parameters.get("ResourceArn")
        tags = self.parameters.get("TagsToAdd", {})
        self.glue_backend.tag_resource(resource_arn, tags)  # type: ignore[arg-type]
        return 201, {}, "{}"

    def untag_resource(self) -> TYPE_RESPONSE:
        resource_arn = self._get_param("ResourceArn")
        tag_keys = self.parameters.get("TagsToRemove")
        self.glue_backend.untag_resource(resource_arn, tag_keys)  # type: ignore[arg-type]
        return 200, {}, "{}"

    def filter_jobs_by_tags(
        self, jobs: List[FakeJob], tags: Dict[str, str]
    ) -> List[str]:
        if not tags:
            return [job.get_name() for job in jobs]
        return [job.get_name() for job in jobs if self.is_tags_match(job.arn, tags)]

    def filter_triggers_by_tags(
        self, triggers: List[FakeTrigger], tags: Dict[str, str]
    ) -> List[str]:
        if not tags:
            return [trigger.get_name() for trigger in triggers]
        return [
            trigger.get_name()
            for trigger in triggers
            if self.is_tags_match(trigger.arn, tags)
        ]

    def is_tags_match(self, resource_arn: str, tags: Dict[str, str]) -> bool:
        glue_resource_tags = self.glue_backend.get_tags(resource_arn)
        mutual_keys = set(glue_resource_tags).intersection(tags)
        for key in mutual_keys:
            if glue_resource_tags[key] == tags[key]:
                return True
        return False

    def create_registry(self) -> str:
        registry_name = self._get_param("RegistryName")
        description = self._get_param("Description")
        tags = self._get_param("Tags")
        registry = self.glue_backend.create_registry(registry_name, description, tags)
        return json.dumps(registry)

    def delete_registry(self) -> str:
        registry_id = self._get_param("RegistryId")
        registry = self.glue_backend.delete_registry(registry_id)
        return json.dumps(registry)

    def get_registry(self) -> str:
        registry_id = self._get_param("RegistryId")
        registry = self.glue_backend.get_registry(registry_id)
        return json.dumps(registry)

    def list_registries(self) -> str:
        registries = self.glue_backend.list_registries()
        return json.dumps({"Registries": registries})

    def create_schema(self) -> str:
        registry_id = self._get_param("RegistryId")
        schema_name = self._get_param("SchemaName")
        data_format = self._get_param("DataFormat")
        compatibility = self._get_param("Compatibility")
        description = self._get_param("Description")
        tags = self._get_param("Tags")
        schema_definition = self._get_param("SchemaDefinition")
        schema = self.glue_backend.create_schema(
            registry_id,
            schema_name,
            data_format,
            compatibility,
            schema_definition,
            description,
            tags,
        )
        return json.dumps(schema)

    def register_schema_version(self) -> str:
        schema_id = self._get_param("SchemaId")
        schema_definition = self._get_param("SchemaDefinition")
        schema_version = self.glue_backend.register_schema_version(
            schema_id, schema_definition
        )
        return json.dumps(schema_version)

    def get_schema_version(self) -> str:
        schema_id = self._get_param("SchemaId")
        schema_version_id = self._get_param("SchemaVersionId")
        schema_version_number = self._get_param("SchemaVersionNumber")

        schema_version = self.glue_backend.get_schema_version(
            schema_id, schema_version_id, schema_version_number
        )
        return json.dumps(schema_version)

    def get_schema_by_definition(self) -> str:
        schema_id = self._get_param("SchemaId")
        schema_definition = self._get_param("SchemaDefinition")
        schema_version = self.glue_backend.get_schema_by_definition(
            schema_id, schema_definition
        )
        return json.dumps(schema_version)

    def put_schema_version_metadata(self) -> str:
        schema_id = self._get_param("SchemaId")
        schema_version_number = self._get_param("SchemaVersionNumber")
        schema_version_id = self._get_param("SchemaVersionId")
        metadata_key_value = self._get_param("MetadataKeyValue")
        schema_version = self.glue_backend.put_schema_version_metadata(
            schema_id, schema_version_number, schema_version_id, metadata_key_value
        )
        return json.dumps(schema_version)

    def get_schema(self) -> str:
        schema_id = self._get_param("SchemaId")
        schema = self.glue_backend.get_schema(schema_id)
        return json.dumps(schema)

    def delete_schema(self) -> str:
        schema_id = self._get_param("SchemaId")
        schema = self.glue_backend.delete_schema(schema_id)
        return json.dumps(schema)

    def update_schema(self) -> str:
        schema_id = self._get_param("SchemaId")
        compatibility = self._get_param("Compatibility")
        description = self._get_param("Description")
        schema = self.glue_backend.update_schema(schema_id, compatibility, description)
        return json.dumps(schema)

    def create_session(self) -> str:
        self.glue_backend.create_session(
            session_id=self.parameters.get("Id"),  # type: ignore[arg-type]
            description=self.parameters.get("Description"),  # type: ignore[arg-type]
            role=self.parameters.get("Role"),  # type: ignore[arg-type]
            command=self.parameters.get("Command"),  # type: ignore[arg-type]
            timeout=self.parameters.get("Timeout"),  # type: ignore[arg-type]
            idle_timeout=self.parameters.get("IdleTimeout"),  # type: ignore[arg-type]
            default_arguments=self.parameters.get("DefaultArguments"),  # type: ignore[arg-type]
            connections=self.parameters.get("Connections"),  # type: ignore[arg-type]
            max_capacity=self.parameters.get("MaxCapacity"),  # type: ignore[arg-type]
            number_of_workers=self.parameters.get("NumberOfWorkers"),  # type: ignore[arg-type]
            worker_type=self.parameters.get("WorkerType"),  # type: ignore[arg-type]
            security_configuration=self.parameters.get("SecurityConfiguration"),  # type: ignore[arg-type]
            glue_version=self.parameters.get("GlueVersion"),  # type: ignore[arg-type]
            tags=self.parameters.get("Tags"),  # type: ignore[arg-type]
            request_origin=self.parameters.get("RequestOrigin"),  # type: ignore[arg-type]
        )
        return ""

    def get_session(self) -> str:
        session_id = self.parameters.get("Id")
        session = self.glue_backend.get_session(session_id)  # type: ignore[arg-type]
        return json.dumps({"Session": session.as_dict()})

    def list_sessions(self) -> str:
        next_token = self._get_param("NextToken")
        max_results = self._get_int_param("MaxResults")
        tags = self._get_param("Tags")
        sessions, next_token = self.glue_backend.list_sessions(
            next_token=next_token, max_results=max_results
        )
        filtered_session_ids = self._filter_sessions_by_tags(sessions, tags)

        return json.dumps(
            dict(
                Ids=[session_id for session_id in filtered_session_ids],
                Sessions=[
                    self.glue_backend.get_session(session_id).as_dict()
                    for session_id in filtered_session_ids
                ],
                NextToken=next_token,
            )
        )

    def _filter_sessions_by_tags(
        self, sessions: List[FakeSession], tags: Dict[str, str]
    ) -> List[str]:
        if not tags:
            return [session.get_id() for session in sessions]
        return [
            session.get_id()
            for session in sessions
            if self.is_tags_match(session.arn, tags)
        ]

    def stop_session(self) -> str:
        session_id = self.parameters.get("Id")
        self.glue_backend.stop_session(session_id)  # type: ignore[arg-type]
        return json.dumps({"Id": session_id})

    def delete_session(self) -> str:
        session_id = self.parameters.get("Id")
        self.glue_backend.delete_session(session_id)  # type: ignore[arg-type]
        return json.dumps({"Id": session_id})

    def batch_get_crawlers(self) -> str:
        crawler_names = self._get_param("CrawlerNames")
        crawlers = self.glue_backend.batch_get_crawlers(crawler_names)
        crawlers_not_found = list(
            set(crawler_names) - set(map(lambda crawler: crawler["Name"], crawlers))
        )
        return json.dumps(
            {
                "Crawlers": crawlers,
                "CrawlersNotFound": crawlers_not_found,
            }
        )

    def batch_get_jobs(self) -> str:
        job_names = self._get_param("JobNames")
        jobs = self.glue_backend.batch_get_jobs(job_names)
        jobs_not_found = list(set(job_names) - set(map(lambda job: job["Name"], jobs)))
        return json.dumps(
            {
                "Jobs": jobs,
                "JobsNotFound": jobs_not_found,
            }
        )

    def get_partition_indexes(self) -> str:
        return json.dumps({"PartitionIndexDescriptorList": []})

    def create_trigger(self) -> str:
        name = self._get_param("Name")
        workflow_name = self._get_param("WorkflowName")
        trigger_type = self._get_param("Type")
        schedule = self._get_param("Schedule")
        predicate = self._get_param("Predicate")
        actions = self._get_param("Actions")
        description = self._get_param("Description")
        start_on_creation = self._get_param("StartOnCreation")
        tags = self._get_param("Tags")
        event_batching_condition = self._get_param("EventBatchingCondition")
        self.glue_backend.create_trigger(
            name=name,
            workflow_name=workflow_name,
            trigger_type=trigger_type,
            schedule=schedule,
            predicate=predicate,
            actions=actions,
            description=description,
            start_on_creation=start_on_creation,
            tags=tags,
            event_batching_condition=event_batching_condition,
        )
        return json.dumps({"Name": name})

    def get_trigger(self) -> str:
        name = self.parameters.get("Name")
        trigger = self.glue_backend.get_trigger(name)  # type: ignore[arg-type]
        return json.dumps({"Trigger": trigger.as_dict()})

    def get_triggers(self) -> str:
        next_token = self._get_param("NextToken")
        dependent_job_name = self._get_param("DependentJobName")
        max_results = self._get_int_param("MaxResults")
        triggers, next_token = self.glue_backend.get_triggers(
            next_token=next_token,
            dependent_job_name=dependent_job_name,
            max_results=max_results,
        )
        return json.dumps(
            dict(
                Triggers=[trigger.as_dict() for trigger in triggers],
                NextToken=next_token,
            )
        )

    def list_triggers(self) -> str:
        next_token = self._get_param("NextToken")
        dependent_job_name = self._get_param("DependentJobName")
        max_results = self._get_int_param("MaxResults")
        tags = self._get_param("Tags")
        triggers, next_token = self.glue_backend.list_triggers(
            next_token=next_token,
            dependent_job_name=dependent_job_name,
            max_results=max_results,
        )
        filtered_trigger_names = self.filter_triggers_by_tags(triggers, tags)
        return json.dumps(
            dict(
                TriggerNames=[trigger_name for trigger_name in filtered_trigger_names],
                NextToken=next_token,
            )
        )

    def batch_get_triggers(self) -> str:
        trigger_names = self._get_param("TriggerNames")
        triggers = self.glue_backend.batch_get_triggers(trigger_names)
        triggers_not_found = list(
            set(trigger_names) - set(map(lambda trigger: trigger["Name"], triggers))
        )
        return json.dumps(
            {
                "Triggers": triggers,
                "TriggersNotFound": triggers_not_found,
            }
        )

    def start_trigger(self) -> str:
        name = self.parameters.get("Name")
        self.glue_backend.start_trigger(name)  # type: ignore[arg-type]
        return json.dumps({"Name": name})

    def stop_trigger(self) -> str:
        name = self.parameters.get("Name")
        self.glue_backend.stop_trigger(name)  # type: ignore[arg-type]
        return json.dumps({"Name": name})

    def delete_trigger(self) -> str:
        name = self.parameters.get("Name")
        self.glue_backend.delete_trigger(name)  # type: ignore[arg-type]
        return json.dumps({"Name": name})
